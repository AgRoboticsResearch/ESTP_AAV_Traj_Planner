from typing import List
import numpy as np
import math
from copy import copy
import matplotlib.path as mpltPath

from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely.affinity import rotate, translate
from shapely.ops import unary_union
from scipy.spatial import ConvexHull


class VehicleModel:
    def __init__(
        self,
        max_steer=0.55,
        wheel_base=1.9,
        axle_to_front=2.85,
        axle_to_back=0.5,
        width=1.48,
        head_out=0.542,
        head_side=0.44,
        car_vertices=None,
        aux_vertices_list=None,
        enable_aux=False,
        enable_plt=True,
    ):
        # Define vehicle parameters
        self.MAX_STEER = max_steer
        self.WHEEL_BASE = wheel_base
        self.AXLE_TO_FRONT = axle_to_front
        self.AXLE_TO_BACK = axle_to_back
        self.WIDTH = width
        self.HEAD_OUT = head_out
        self.HEAD_SIDE = head_side
        self.MAX_CURVATURE = math.tan(max_steer) / wheel_base

        # Initialize geometry
        self.car_vertices = car_vertices if car_vertices else []
        self.aux_vertices_list = aux_vertices_list if aux_vertices_list else []
        self.enable_aux = enable_aux
        self.enable_plt = enable_plt

        self.car_poly = []
        self.aux_polys = []
        self.plt_car_poly = []
        self.plt_aux_polys = []

        self.initialize_geometry()

    def initialize_geometry(self):
        if self.car_vertices:
            self.generate_custom_car_polygon(self.car_vertices)
        else:
            self.generate_rectangle_car_polygon()

        self.car_poly = orient(self.car_poly, sign=-1)  # ensure clockwise orientation
        if self.enable_plt:
            self.plt_car_poly = mpltPath.Path(
                self.car_poly.exterior.coords
            )  # for matplot

        if self.enable_aux and self.aux_vertices_list:
            self.generate_custom_aux_polygons(self.aux_vertices_list)

    def generate_custom_car_polygon(self, vertices):
        assert len(vertices) >= 3, "Car polygon must have at least 3 vertices!"
        self.reassign_first_point(vertices)
        self.car_poly = Polygon(vertices)

    def generate_rectangle_car_polygon(self):
        self.car_vertices = [
            [-self.AXLE_TO_BACK, self.WIDTH / 2],  # top left point
            [-self.AXLE_TO_BACK, -self.WIDTH / 2],
            [self.AXLE_TO_FRONT, -self.WIDTH / 2],
            [self.AXLE_TO_FRONT, self.WIDTH / 2],
            [-self.AXLE_TO_BACK, self.WIDTH / 2],
        ]
        self.car_poly = Polygon(self.car_vertices)

    def generate_custom_aux_polygons(self, vertices_list):
        if len(vertices_list) == 0:
            return

        for vertices in vertices_list:
            # if vertices is a tuple of ([x, y], dx, dy)
            if (
                len(vertices) == 3
                and not isinstance(vertices[1], List)
                and not isinstance(vertices[2], List)
            ):
                top_left_coord, dx, dy = vertices[0], vertices[1], vertices[2]
                p1 = np.array(top_left_coord)
                p2 = np.array([top_left_coord[0] + dx, top_left_coord[1]])
                p3 = np.array([top_left_coord[0] + dx, top_left_coord[1] - dy])
                p4 = np.array([top_left_coord[0], top_left_coord[1] - dy])
                poly_points = np.array([p1, p2, p3, p4])

            else:  # else if vertices is a list of [x, y] points
                poly_points = vertices
                self.reassign_first_point(poly_points)

            aux_poly = orient(
                Polygon(poly_points), sign=-1
            )  # ensure clockwise orientation
            self.aux_polys.append(aux_poly)

        if self.enable_plt:
            self.plt_aux_polys = [
                mpltPath.Path(poly.exterior.coords) for poly in self.aux_polys
            ]

    def draw_car(
        self, plt, x, y, yaw, car_color="orange", aux_color="orange", alpha=0.1
    ):
        car_poly_in_odom, aux_polys_in_odom = self.get_car_aux_polys_in_odom(x, y, yaw)

        car_x, car_y = car_poly_in_odom.exterior.xy
        plt.plot(car_x, car_y, color="k", alpha=alpha)
        plt.fill(car_x, car_y, color=car_color, alpha=alpha)

        if self.enable_aux and len(aux_polys_in_odom) > 0:
            for aux_poly in aux_polys_in_odom:
                aux_x, aux_y = aux_poly.exterior.xy
                plt.plot(aux_x, aux_y, color="k", alpha=alpha)
                plt.fill(aux_x, aux_y, color=aux_color, alpha=alpha)

    def draw_aux(self, plt, x, y, yaw, aux_color="orange", alpha=0.1):
        _, aux_polys_in_odom = self.get_car_aux_polys_in_odom(x, y, yaw)

        if self.enable_aux and len(aux_polys_in_odom) > 0:
            for aux_poly in aux_polys_in_odom:
                aux_x, aux_y = aux_poly.exterior.xy
                plt.plot(aux_x, aux_y, color="k", alpha=alpha)
                plt.fill(aux_x, aux_y, color=aux_color, alpha=alpha)

    def draw_conservative_convex_poly(self, plt, x, y, yaw, color="orange", alpha=0.1):
        convex_poly_in_odom = self.get_conservative_convex_poly_in_odom(x, y, yaw)
        convex_hull_x, convex_hull_y = convex_poly_in_odom.exterior.xy
        plt.plot(convex_hull_x, convex_hull_y, color="k", alpha=alpha)
        plt.fill(convex_hull_x, convex_hull_y, color=color, alpha=alpha)

    def get_car_poly_in_body(self, buffer=0.0):
        car_poly = copy(self.car_poly)
        car_poly = car_poly.buffer(buffer)

        return car_poly

    def get_aux_polys_in_body(self, buffer=0.0):
        aux_polys = copy(self.aux_polys)
        for aux_poly in aux_polys:
            aux_poly = aux_poly.buffer(buffer)

        return aux_polys

    def get_car_aux_polys_in_odom(self, x, y, yaw_in_rad, buffer=0.0):
        car_poly_in_body = self.get_car_poly_in_body(buffer)
        aux_polys_in_body = self.get_aux_polys_in_body(buffer)

        car_poly_in_odom = self.transform_polygon(car_poly_in_body, x, y, yaw_in_rad)

        aux_polys_in_odom = []
        if self.enable_aux and len(aux_polys_in_body) > 0:
            for aux_poly in aux_polys_in_body:
                aux_poly_in_odom = self.transform_polygon(aux_poly, x, y, yaw_in_rad)
                aux_polys_in_odom.append(aux_poly_in_odom)

        return car_poly_in_odom, aux_polys_in_odom

    def get_conservative_convex_poly_in_odom(self, x, y, yaw_in_rad, buffer=0.0):
        convex_poly_in_body = self.get_conservative_convex_poly_in_body(buffer)
        convex_poly_in_odom = self.transform_polygon(
            convex_poly_in_body, x, y, yaw_in_rad
        )

        return convex_poly_in_odom

    def get_path_polys(self, path: np.ndarray, skip=1):
        car_path_polys = [
            self.transform_polygon(self.car_poly, *pose) for pose in path[::skip]
        ]
        car_path_polys = unary_union(car_path_polys)

        aux_path_polys = []
        if self.enable_aux and len(self.aux_polys) > 0:
            for aux_poly in self.aux_polys:
                aux_path_polys += [
                    self.transform_polygon(aux_poly, *pose) for pose in path[::skip]
                ]
            aux_path_polys = unary_union(aux_path_polys)

        return car_path_polys, aux_path_polys

    def get_car_vertices_in_body(self, buffer=0.0):
        car_vertices_in_body = np.asarray(
            self.get_car_poly_in_body(buffer).exterior.xy
        ).T

        return car_vertices_in_body

    def get_aux_vertices_in_body(self, buffer=0.0):
        aux_vertices_list_in_body = []
        for aux_poly in self.get_aux_polys_in_body(buffer):
            coords_in_body = np.asarray(aux_poly.exterior.xy).T
            aux_vertices_list_in_body.append(coords_in_body)

        return aux_vertices_list_in_body

    def get_car_convex_poly_in_body(self, buffer=0.0):
        car_vertices = self.get_car_vertices_in_body()
        convex_hull = ConvexHull(car_vertices)
        convex_vertices = convex_hull.points[convex_hull.vertices].tolist()
        self.reassign_first_point(convex_vertices)
        convex_poly = Polygon(convex_vertices)
        convex_poly = orient(convex_poly, sign=-1)  # ensure clockwise orientation
        convex_poly = convex_poly.buffer(buffer)

        return convex_poly

    def get_aux_convex_poly_in_body(self, buffer=0.0):
        aux_vertices_list = self.get_aux_vertices_in_body()
        convex_poly_list = []
        for aux_vertices in aux_vertices_list:
            convex_hull = ConvexHull(aux_vertices)
            convex_vertices = convex_hull.points[convex_hull.vertices].tolist()
            self.reassign_first_point(convex_vertices)
            convex_poly = Polygon(convex_vertices)
            convex_poly = orient(convex_poly, sign=-1)  # ensure clockwise orientation
            convex_poly = convex_poly.buffer(buffer)
            convex_poly_list.append(convex_poly)

        return convex_poly_list

    def get_conservative_convex_poly_in_body(self, buffer=0.0):
        car_vertices = self.get_car_vertices_in_body()
        aux_vertices_list = self.get_aux_vertices_in_body()
        convex_hull = ConvexHull(np.vstack([car_vertices, *aux_vertices_list]))
        convex_poly = Polygon(convex_hull.points[convex_hull.vertices])
        convex_poly = orient(convex_poly, sign=-1)  # ensure clockwise orientation
        convex_poly = convex_poly.buffer(buffer)

        return convex_poly

    def get_car_convex_poly_vertices_in_body(self, buffer=0.0):
        convex_poly = self.get_car_convex_poly_in_body(buffer)
        convex_vertices = np.asarray(convex_poly.exterior.xy).T

        return convex_vertices

    def get_aux_convex_poly_vertices_list_in_body(self, buffer=0.0):
        convex_poly_list = self.get_aux_convex_poly_in_body(buffer)
        convex_vertices_list = [
            np.asarray(convex_poly.exterior.xy).T for convex_poly in convex_poly_list
        ]

        return convex_vertices_list

    def get_conservative_convex_poly_vertices_in_body(self, buffer=0.0):
        convex_poly = self.get_conservative_convex_poly_in_body(buffer)
        convex_vertices = np.asarray(convex_poly.exterior.xy).T

        return convex_vertices

    def get_minimum_turn_radius(self, max_steer_angle=None):
        if max_steer_angle is None:
            return 1.0 / self.MAX_CURVATURE
        else:
            return self.WHEEL_BASE / math.tan(max_steer_angle)

    def get_steer_angle(self, curvature):
        steer_angle = np.clip(
            math.atan(self.WHEEL_BASE * curvature), -self.MAX_STEER, self.MAX_STEER
        )

        return steer_angle

    @staticmethod
    def transform_polygon(poly: Polygon, x, y, yaw):
        new_poly = rotate(poly, yaw, use_radians=True, origin=(0, 0))
        new_poly = translate(new_poly, x, y)

        return new_poly

    @staticmethod
    def reassign_first_point(vertices: List):
        if vertices[0] == vertices[-1]:
            vertices.pop()

        # use top left as the first point
        first_point_idx = min(
            range(len(vertices)), key=lambda i: (vertices[i][0], -vertices[i][1])
        )
        vertices[:] = vertices[first_point_idx:] + vertices[:first_point_idx]

        if vertices[0] != vertices[-1]:
            vertices.append(vertices[0])