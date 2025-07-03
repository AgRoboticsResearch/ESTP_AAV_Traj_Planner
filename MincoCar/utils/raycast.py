import numpy as np


def mod(value: float, modulus: float) -> float:
    return np.fmod(np.fmod(value, modulus) + modulus, modulus)


def int_bound(s: float, ds: float) -> float:
    """Find the smallest positive t such that s+t*ds is an integer"""
    if ds < 0:
        return int_bound(-s, -ds)
    else:
        s = mod(s, 1)
        # problem is now s+t*ds = 1
        return (1 - s) / ds


class RayCaster2D:
    """
    RayCaster2D Class

    From "A Fast Voxel Traversal Algorithm for Ray Tracing"
    by John Amanatides and Andrew Woo, 1987
    <http://www.cse.yorku.ca/~amana/research/grid.pdf>
    <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.3443>
    Extensions to the described algorithm:
        • Imposed a distance limit.
        • The face passed through to reach the current cube is provided to
        the callback.

    The foundation of this algorithm is a parameterized representation of
    the provided ray,
                    origin + t * direction,
    except that t is not actually stored; rather, at any given point in the
    traversal, we keep track of the *greater* t values which we would have
    if we took a step sufficient to cross a cube boundary along that axis
    (i.e. change the integer part of the coordinate) in the variables
    tMaxX, tMaxY, and tMaxZ.

    """

    def __init__(self):
        pass

    def set_input(self, start: np.ndarray, end: np.ndarray) -> bool:
        self.start = start
        self.end = end

        # voxel containing origin point
        self.x, self.y = np.floor(self.start).astype(int)
        self.endX, self.endY = np.floor(self.end).astype(int)
        direction = end - start

        # Break out direction vector
        # dx = self.endX - self.x
        # dy = self.endY - self.y
        self.dx = direction[0]
        self.dy = direction[1]

        # Direcion to increment x,y when stepping
        self.stepX = np.sign(self.dx)
        self.stepY = np.sign(self.dy)

        # See description above. The initial values depend on the fractional
        # part of the origin
        self.tMaxX = int_bound(self.start[0], self.dx) if self.dx != 0 else float("inf")
        self.tMaxY = int_bound(self.start[1], self.dy) if self.dy != 0 else float("inf")

        # The change in t when taking a step (always positive)
        self.tDeltaX = float(self.stepX) / self.dx if self.dx != 0 else float("inf")
        self.tDeltaY = float(self.stepY) / self.dy if self.dy != 0 else float("inf")

        # Avoid an infinte loop
        if self.stepX == 0 and self.stepY == 0:
            return False
        return True

    def step(self, ray_pt: np.ndarray) -> bool:
        ray_pt[:] = np.array([self.x, self.y], dtype=float)

        if self.x == self.endX and self.y == self.endY:
            return False

        if self.tMaxX < self.tMaxY:
            self.x += self.stepX
            self.tMaxX += self.tDeltaX
        else:
            self.y += self.stepY
            self.tMaxY += self.tDeltaY

        return True
