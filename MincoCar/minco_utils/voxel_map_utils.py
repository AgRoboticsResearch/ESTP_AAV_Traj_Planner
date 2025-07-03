from typing import List, Tuple, Union
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, box, mapping
import shapely.affinity
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy
import rasterio
from rasterio.features import rasterize

from .raycast import RayCaster2D
from .polygon_utils import invert_polygon

from pymincocar import VoxelMap

def create_3d_voxel_map_from_2d_occupancy(
    occupancy_map_2d: np.ndarray,
    origin: np.ndarray = np.zeros(3),
    voxel_scale: float = 1.0,
):
    height, width = occupancy_map_2d.shape
    z_dim = 1
    map_size = np.array([width, height, z_dim], dtype=np.int32)
    map_origin = np.array([origin[0], origin[1], 0.0])

    voxel_map = VoxelMap(map_size, map_origin, voxel_scale)

    # set occupancy
    occupied_positions = np.where(occupancy_map_2d == 1)
    for y, x in zip(*occupied_positions):
        for z in range(z_dim):
            voxel_map.setOccupied(np.array([x, y, z], dtype=np.int32))

    return voxel_map


def create_3d_voxel_map_from_grid_map(grid_map: "GridMap2D"):
    occupancy_map_2d = grid_map.get_grid().T
    origin = grid_map.get_origin()
    voxel_scale = grid_map.get_resolution()

    voxel_map = create_3d_voxel_map_from_2d_occupancy(
        occupancy_map_2d, origin, voxel_scale
    )

    return voxel_map


def draw_3d_voxel_map_in_2d(ax: plt.Axes, voxel_map: VoxelMap):
    map_voxels = voxel_map.getVoxels()
    map_size = voxel_map.getGridSize()
    map_origin = voxel_map.getOrigin()
    map_voxel_scale = voxel_map.getScale()

    map_2d = np.zeros((map_size[1], map_size[0]), dtype=int)
    for idx, state in enumerate(map_voxels):
        if state == 1:
            x = idx % map_size[0]
            y = (idx // map_size[0]) % map_size[1]
            map_2d[y, x] = 1

    extent = [
        map_origin[0],
        map_origin[0] + map_size[0] * map_voxel_scale,
        map_origin[1],
        map_origin[1] + map_size[1] * map_voxel_scale,
    ]

    ax.imshow(map_2d, cmap="Greys", origin="lower", extent=extent)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

def convert_voxel_map_to_2d(voxel_map):
    # Get voxel data and map parameters
    voxels_2d = voxel_map.getVoxels()
    grid_size = voxel_map.getGridSize()
    origin = voxel_map.getOrigin()
    resolution = voxel_map.getResolution()

    # Create a 2D array to store the map data
    map_2d = np.zeros((grid_size[1], grid_size[0]), dtype=int)
    
    # Iterate over the voxel data and fill in the 2D map
    for idx, state in enumerate(voxels_2d):
        x = idx % grid_size[0]
        y = (idx // grid_size[0]) % grid_size[1]
        # Invert the state: if occupied (state > 0), set to 0 (free), otherwise set to 1 (occupied)
        map_2d[y, x] = 0 if state > 0 else 1
    return map_2d, grid_size, origin, resolution

def create_voxel_map_from_2d_occupancy(
    occupancy_map_2d: np.ndarray,
    origin: np.ndarray = np.zeros(3),
    resolution: float = 1.0,
    dilate: int = 1,
):
    height, width = occupancy_map_2d.shape
    z_dim = 1
    grid_size = np.array([width, height, z_dim], dtype=np.int32)
    map_origin = np.array([origin[0], origin[1], 0.0])

    voxel_map = VoxelMap(grid_size, map_origin, resolution)

    # set 2d occupancy
    occupied_positions = np.where(occupancy_map_2d == 1)
    for y, x in zip(*occupied_positions):
        res = voxel_map.setOccupancy(np.array([x, y], dtype=np.int32), 1)
        if res is False:
            raise ValueError("Error in setting occupancy.")
    
    # dilate
    voxel_map.dilate(int(dilate))
            
    return voxel_map


def create_voxel_map_from_grid_map(grid_map: "GridMap2D", dilate: int = 1):
    occupancy_map_2d = grid_map.get_grid().T
    origin = grid_map.get_origin()
    voxel_scale = grid_map.get_resolution()

    voxel_map = create_voxel_map_from_2d_occupancy(
        occupancy_map_2d, origin, voxel_scale, dilate
    )

    return voxel_map


def draw_voxel_map_in_2d(ax: plt.Axes, voxel_map: VoxelMap):
    voxels_2d = voxel_map.getVoxels()
    grid_size = voxel_map.getGridSize()
    origin = voxel_map.getOrigin()
    resolution = voxel_map.getResolution()

    map_2d = np.zeros((grid_size[1], grid_size[0]), dtype=int)
    for idx, state in enumerate(voxels_2d):
        if state > 0:
            x = idx % grid_size[0]
            y = (idx // grid_size[0]) % grid_size[1]
            map_2d[y, x] = state
            
    extent = [
        origin[0],
        origin[0] + grid_size[0] * resolution,
        origin[1],
        origin[1] + grid_size[1] * resolution,
    ]

    cmap = plt.cm.get_cmap("Greys", 3)
    ax.imshow(map_2d, cmap=cmap, alpha=0.7, origin="lower", extent=extent)
    ax.set_xlabel("x")
    ax.set_ylabel("y")    
class FeatureMap2D:
    """
    Feature Map in 2D Class
    """

    def __init__(
        self,
        x_min: float = 0.0,
        y_min: float = 0.0,
        x_max: float = 20.0,
        y_max: float = 20.0,
    ):
        # create boundary
        self.create_boundary(x_min, y_min, x_max, y_max)

        # obstacle polygons
        self.obs_polygons = []

    def create_boundary(self, x_min: float, y_min: float, x_max: float, y_max: float):
        if x_min >= x_max:
            raise ValueError("x_min should be smaller than x_max")
        if y_min >= y_max:
            raise ValueError("y_min should be smaller than y_max")

        self.x_min, self.y_min = float(x_min), float(y_min)
        self.x_max, self.y_max = float(x_max), float(y_max)
        self.map_size_x = self.x_max - self.x_min
        self.map_size_y = self.y_max - self.y_min
        self.boundary_box = box(self.x_min, self.y_min, self.x_max, self.y_max)

    def cast_polygon_within_box(self, polygon: Polygon) -> Polygon:
        """Cast the polygon within the boundary box"""
        cast_polygon = polygon.intersection(self.boundary_box)
        return cast_polygon

    def add_polygon(self, polygon: Polygon, is_relative: bool = False):
        if not isinstance(polygon, Polygon):
            raise TypeError(f"{polygon} is not a Polygon")

        if is_relative:
            translated_polygon = shapely.affinity.translate(
                polygon, self.x_min, self.y_min
            )
            casted_polygon = self.cast_polygon_within_box(translated_polygon)
        else:
            casted_polygon = self.cast_polygon_within_box(polygon)

        if not casted_polygon.is_empty:
            self.obs_polygons.append(casted_polygon)

    def add_polygons(self, polygon_list: List[Polygon], is_relative: bool = False):
        for polygon in polygon_list:
            self.add_polygon(polygon, is_relative)

    def pop(self):
        self.obs_polygons.pop()

    def reset(self):
        self.obs_polygons = []

    def get_polys(self) -> List[Polygon]:
        return self.obs_polygons

    def get_polys_num(self) -> int:
        return len(self.obs_polygons)

    def get_boundary(self) -> Tuple:
        return (self.x_min, self.y_min, self.x_max, self.y_max)

    def draw(
        self,
        ax: plt.Axes,
        fill: bool = True,
        alpha: float = 1.0,
        facecolor: str = "k",
        edgecolor: str = "k",
        hatch: str = None,
        linewidth: float = 1.5,
        title: str = "2D Feature Map",
    ):
        # plot boundary
        ax.plot(*self.boundary_box.exterior.xy, color="k", linestyle="-", linewidth=2.0)
        # patch = patches.Polygon(list(self.boundary_box.exterior.coords), color='gray', fill=True, alpha=0.1)
        # ax.add_patch(patch)

        # plot obstacle polygons
        for polygon in self.obs_polygons:
            patch = patches.Polygon(
                list(polygon.exterior.coords),
                fill=fill,
                alpha=alpha,
                facecolor=facecolor,
                edgecolor=edgecolor,
                hatch=hatch,
                linewidth=linewidth,
            )
            ax.add_patch(patch)

        # leave a margin for better visualization
        margin_x = self.map_size_x * 0.05
        margin_y = self.map_size_y * 0.05
        ax.set_xlim([self.x_min - margin_x, self.x_max + margin_x])
        ax.set_ylim([self.y_min - margin_y, self.y_max + margin_y])
        ax.set_aspect(1)
        ax.set_title(title)
        # ax.set_xlabel("x")
        # ax.set_ylabel("y")

    def to_gridmap(
        self, resolution: float = 0.2, use_raster: bool = True
    ) -> "GridMap2D":
        if resolution <= 0:
            raise ValueError("Resolution should be positive.")

        grid_map = GridMap2D(
            self.x_min,
            self.y_min,
            self.map_size_x,
            self.map_size_y,
            resolution,
        )

        # Initialize an empty grid
        grid_size = grid_map.get_grid_size()
        grid = np.zeros(grid_size, dtype=int)

        # 1) Use rasterio (fast but inaccurate)
        if use_raster:
            # Create a transform
            transform = rasterio.transform.from_bounds(
                self.y_min,
                self.x_max,
                self.y_max,
                self.x_min,
                grid_size[1],
                grid_size[0],
            )

            # Rasterize each polygon
            for polygon in self.obs_polygons:
                # The rasterize function converts the polygon into a binary mask
                mask = rasterio.features.rasterize(
                    [(mapping(invert_polygon(polygon)), 1)],
                    out_shape=grid_size,
                    transform=transform,
                    fill=0,
                    all_touched=True,
                    dtype=np.uint8,
                )
                grid += mask

        # 2) manual (slow but accurate)
        else:
            multi_polys = MultiPolygon(self.obs_polygons)

            for x in np.arange(
                self.x_min + resolution / 2,
                self.x_max - resolution / 2 + 1e-5,
                resolution,
            ):
                for y in np.arange(
                    self.y_min + resolution / 2,
                    self.y_max - resolution / 2 + 1e-5,
                    resolution,
                ):
                    if multi_polys.contains(Point(x, y)):
                        x_id = np.floor((x - self.x_min) / resolution).astype(int)
                        y_id = np.floor((y - self.y_min) / resolution).astype(int)
                        grid[x_id][y_id] += 1

        # grid is a grid map where each cell contains 1 if it is inside a polygon and 0 otherwise
        np.putmask(grid, grid > 0, 1)

        grid_map.set_value_from_grid(grid)
        return grid_map


class GridMap2D:
    """
    Grid Map in 2D Class
    """

    def __init__(
        self,
        origin_x: float = 0.0,  # bottom-left
        origin_y: float = 0.0,
        map_size_x: float = 20.0,
        map_size_y: float = 20.0,
        resolution: float = 0.2,
    ):
        if map_size_x <= 0 or map_size_y <= 0:
            raise ValueError("Map size should be positive.")

        if resolution <= 0:
            raise ValueError("Resolution should be positive.")

        self.origin = np.array([origin_x, origin_y], dtype=float)
        self.resolution = resolution

        self.grid_size = (
            np.ceil(map_size_x / self.resolution).astype(int),
            np.ceil(map_size_y / self.resolution).astype(int),
        )

        if self.grid_size[0] <= 0 or self.grid_size[1] <= 0:
            raise ValueError("Grid size error.")

        self.map_size = (
            float(self.grid_size[0] * resolution),
            float(self.grid_size[1] * resolution),
        )

        # initialize grid map
        self.grid = np.zeros(self.grid_size, dtype=float)

    def __getitem__(self, indices: List) -> Union[float, None]:
        if self.grid is None:
            return

        if len(indices) != 2:
            raise ValueError("Indices should have size 2.")

        if indices[0] > self.grid.shape[0] - 1 or indices[1] > self.grid.shape[1] - 1:
            raise IndexError(f"{indices} out of index range of grid.")

        return self.grid[indices[0]][indices[1]]

    def set_value_from_grid(self, grid: np.ndarray) -> bool:
        if grid.shape != self.grid.shape:
            print(
                f"Error! Input has a dimension of {grid.shape} while current grid has a dimension of {self.grid.shape}."
            )
            return False

        self.grid = copy.copy(grid).astype(float)
        return True

    def set_value_from_pos(self, pos: np.ndarray, val: float) -> bool:
        index = self.pos_to_index(pos)
        if not index:
            return False

        flag = self.set_value_from_index(index, val)
        return flag

    def set_value_from_index(self, index: Tuple, val: float) -> bool:
        if not self.is_index_in_map(index):
            return False

        self.grid[index[0]][index[1]] = val
        return True

    def get_resolution(self) -> float:
        return self.resolution

    def get_origin(self) -> np.ndarray:
        return self.origin

    def get_map_size(self) -> Tuple:
        return self.map_size

    def get_grid_size(self) -> Tuple:
        return self.grid.shape

    def get_grid(self) -> np.ndarray:
        return self.grid

    def pos_to_index(self, pos: np.ndarray) -> Union[Tuple, None]:
        index = (
            np.floor((pos[0] - self.origin[0]) / self.resolution).astype(int),
            np.floor((pos[1] - self.origin[1]) / self.resolution).astype(int),
        )
        if self.is_index_in_map(index):
            return index
        return

    def index_to_pos(self, index: Tuple) -> np.ndarray:
        if not self.is_index_in_map(index):
            return

        pos = np.array([self.origin[0], self.origin[1]])
        pos[0] += (index[0] + 0.5) * self.resolution
        pos[1] += (index[1] + 0.5) * self.resolution
        return pos

    def is_index_in_map(self, index: Tuple) -> bool:
        result = (
            index[0] >= 0
            and index[0] <= self.grid_size[0] - 1
            and index[1] >= 0
            and index[1] <= self.grid_size[1] - 1
        )
        return result

    def is_pos_in_map(self, pos: np.ndarray) -> bool:
        result = (
            pos[0] >= self.origin[0]
            and pos[0] <= self.origin[0] + self.map_size[0]
            and pos[1] >= self.origin[1]
            and pos[1] <= self.origin[1] + self.map_size[1]
        )
        return result

    def is_state_valid(self, pos: np.ndarray) -> bool:
        index = self.pos_to_index(pos)
        if index:
            return self.grid[index[0]][index[1]] == 0
        return False

    def is_segment_valid(
        self, p0: np.ndarray, p1: np.ndarray, max_dist: float = float("inf")
    ) -> bool:
        if (not self.is_state_valid(p0)) or (not self.is_state_valid(p1)):
            return False

        diff = p1 - p0
        dist = np.linalg.norm(diff)
        if dist > max_dist:
            return False

        raycaster = RayCaster2D()
        need_ray = raycaster.set_input(
            p0 / self.resolution, p1 / self.resolution
        )  # (ray start, ray end)
        if not need_ray:
            return True

        half = np.array([0.5, 0.5], dtype=float)
        ray_pt = np.zeros(2, dtype=float)
        if not raycaster.step(ray_pt):  # skip the ray start point
            return True

        while raycaster.step(ray_pt):
            current_pt = (ray_pt + half) * self.resolution
            if not self.is_state_valid(current_pt):
                return False

        return True

    def draw(
        self,
        ax: plt.Axes,
        use_index: bool = True,
        cmap: str = "gray_r",
        alpha: float = 0.85,
        title: str = "2D Grid Map",
    ):
        if use_index:
            ax.imshow(self.grid.T, cmap=cmap, origin="lower", alpha=alpha)
            ax.set_xlabel("x id")
            ax.set_ylabel("y id")
        else:
            boundary_box = box(
                self.origin[0],
                self.origin[1],
                self.origin[0] + self.map_size[0],
                self.origin[1] + self.map_size[1],
            )
            ax.plot(*boundary_box.exterior.xy, color="k", linestyle="-", linewidth=1.5)

            x = np.array(
                [
                    self.origin[0] + i * self.resolution
                    for i in range(self.grid.shape[0] + 1)
                ]
            )
            y = np.array(
                [
                    self.origin[1] + j * self.resolution
                    for j in range(self.grid.shape[1] + 1)
                ]
            )
            ax.pcolormesh(
                x,
                y,
                self.grid.T,
                cmap=cmap,
                edgecolor="k",
                linewidth=0.01,
                alpha=alpha,
            )

            # leave a margin for better visualization
            margin_x = self.map_size[0] * 0.05
            margin_y = self.map_size[1] * 0.05
            ax.set_xlim(
                [
                    self.origin[0] - margin_x,
                    self.origin[0] + self.map_size[0] + margin_x,
                ]
            )
            ax.set_ylim(
                [
                    self.origin[1] - margin_y,
                    self.origin[1] + self.map_size[1] + margin_y,
                ]
            )
            ax.set_xlabel("x")
            ax.set_ylabel("y")

        ax.set_aspect(1)
        ax.set_title(title)
        ax.invert_yaxis()
