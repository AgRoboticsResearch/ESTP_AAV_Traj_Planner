import numpy as np
import shapely.affinity
from shapely.geometry import Polygon, Point, box
from shapely.ops import unary_union


def create_box(width=1.0, height=1.0, center=[0.0, 0.0], angle=0.0, is_deg=True):
    """Create a box shape"""
    if width <= 0 or height <= 0:
        raise ValueError("Both width and height should be positive.")

    rectangle = box(-0.5, -0.5, 0.5, 0.5)
    rectangle_scale = shapely.affinity.scale(rectangle, width, height)
    if not is_deg:
        angle = np.rad2deg(angle)
    rectangle_scale_rot = shapely.affinity.rotate(rectangle_scale, angle)
    rectangle_scale_rot_trans = shapely.affinity.translate(
        rectangle_scale_rot, center[0], center[1]
    )
    
    return rectangle_scale_rot_trans


def create_circle(center=[0, 0], radius=1):
    """Create an circle shape"""
    if radius < 0:
        raise ValueError("The radius should be positive.")
    
    return Point(center[:2]).buffer(radius)


def create_ellipse(a=1, b=1, center=[0, 0], angle=0, is_deg=True):
    """Create an ellipse shape"""
    if a <= 0 or b <= 0:
        raise ValueError("Both a and b should be positive.")

    circle = create_circle([0, 0], 1)
    circle_scale = shapely.affinity.scale(circle, a, b)
    if not is_deg:
        angle = np.rad2deg(angle)
    circle_scale_rot = shapely.affinity.rotate(circle_scale, angle)
    circle_scale_rot_trans = shapely.affinity.translate(
        circle_scale_rot, center[0], center[1]
    )
    return circle_scale_rot_trans


def create_polygon(vertices):
    """Create a polygon shape"""
    if len(vertices) < 3:
        raise ValueError("The minimum vertices number is 3.")

    try:
        polygon = Polygon(vertices)
        if not polygon.is_valid:
            raise ValueError("The provided vertices do not form a valid polygon.")
        return polygon
    except Exception as e:
        raise Exception(f"Cannot create a polygon with the given vertices! Error: {e}")


def invert_polygon(polygon):
    """Invert a polygon"""
    inverted_vertices = [(y, x) for x, y in polygon.exterior.coords]
    inverted_polygon = Polygon(inverted_vertices)

    return inverted_polygon


def inflate_polygon(polygon, inflation=1.0):
    """Inflate a polygon"""
    if inflation < 0:
        raise ValueError("Inflation should be positive.")
    
    return polygon.buffer(inflation)


def combine_polygons(polygon_list):
    """Combine polygons using convex hull"""
    if len(polygon_list) == 0:
        return Polygon()

    unioned_polygon = unary_union(polygon_list)
    convex_hull = unioned_polygon.convex_hull

    return convex_hull
