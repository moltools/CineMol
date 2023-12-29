"""
This module contains classes for representing 3D shapes.
"""
from dataclasses import dataclass
from enum import Enum, auto
import math
import typing as ty

from cinemol.geometry import Vector3D, Point3D

import numpy as np # TODO: Remove dependency on numpy.

# ==============================================================================
# Helper functions
# ==============================================================================

def gram_schmidt(n: Vector3D) -> ty.Tuple[Vector3D, Vector3D]:
    """
    Generate two vectors that are orthogonal to the given vector.
    
    :param Vector3D n: The vector to generate orthogonal vectors for.
    :return: Two orthogonal vectors to the given vector.
    :rtype: ty.Tuple[Vector3D, Vector3D]
    """
    v = Vector3D.create_random()
    v = v.subtract(n.multiply(v.dot(n))).normalize()
    w = n.cross(v)
    return v, w

# ==============================================================================
# Shapes
# ==============================================================================

@dataclass 
class Sphere:
    """
    A sphere.

    :param np.ndarray center: The center of the sphere, with shape (3,).
    :param float radius: The radius of the sphere.
    """
    center: np.ndarray
    radius: float

    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the sphere.
        
        :param int res: The resolution of the sphere.
        :return: The points on the surface of the sphere, with shape (N, 3).
        :rtype: np.ndarray
        """
        phi = np.linspace(0, 2 * math.pi, res)
        theta = np.linspace(0, math.pi, res)
        phi, theta = np.meshgrid(phi, theta)
        x = self.center[0] + self.radius * np.sin(theta) * np.cos(phi)
        y = self.center[1] + self.radius * np.sin(theta) * np.sin(phi)
        z = self.center[2] + self.radius * np.cos(theta)
        m = np.array([x.flatten(), y.flatten(), z.flatten()])
        return m.T
    
    def point_is_inside(self, point: np.ndarray) -> bool:
        """
        Check if a point is inside the sphere.
        
        :param np.ndarray point: The point to check, with shape (3,).
        :return: True if the point is inside the sphere, False otherwise.
        :rtype: bool
        """
        return np.linalg.norm(point - self.center) <= self.radius
    
@dataclass
class Circle3D:
    """
    A circle in 3D.

    :param np.ndarray center: The center of the circle, with shape (3,).
    :param float radius: The radius of the circle.
    :param np.ndarray normal: The normal vector of the circle, with shape (3,).
    """
    center: np.ndarray
    radius: float
    normal: np.ndarray

    def generate_points_on_circumference(self, res: int) -> np.ndarray:
        """
        Generate points on the circumference of the circle.
        
        :param int res: The resolution of the circle.
        :return: The points on the circumference of the circle, with shape (N, 3).
        :rtype: np.ndarray
        """
        normal = Vector3D(*self.normal)
        normal = normal.normalize()
        v, w = gram_schmidt(normal)
        v = np.array([v.x, v.y, v.z])
        w = np.array([w.x, w.y, w.z])
        angles = np.linspace(0, 2 * np.pi, res)
        return self.center + self.radius * np.outer(np.cos(angles), v) + self.radius * np.outer(np.sin(angles), w)
    
    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the circle.
        
        :param int res: The resolution of the circle.
        :return: The points on the surface of the circle, with shape (N, 3).
        :rtype: np.ndarray
        """
        points = []

        # For each radius, generate points on the circumference.
        radii = np.linspace(0, self.radius, res)
        for radius in radii:
            cap = Circle3D(self.center, radius, self.normal)
            points.extend(cap.generate_points_on_circumference(res))

        return points

@dataclass 
class Line3D:
    """
    A line in 3D.
    
    :param np.ndarray start: The start point of the line, with shape (3,).
    :param np.ndarray end: The end point of the line, with shape (3,).
    """
    start: np.ndarray
    end: np.ndarray
    
    def generate_points_along_line(self, res: int) -> np.ndarray:
        """
        Generate points along the line.

        :param int res: The resolution of the line.
        :return: The points along the line, with shape (N, 3).
        :rtype: np.ndarray
        """
        return np.linspace(self.start, self.end, res)
    
    def distance_to_line(self, point: np.ndarray) -> float:
        """
        Compute the distance from a point to the line.

        :param np.ndarray point: The point, with shape (3,).
        :return: The distance from the point to the line.
        :rtype: float
        """
        # Normalized tangent vector.
        d = np.divide(self.end - self.start, np.linalg.norm(self.end - self.start))

        # Signed parallel distance components.
        s = np.dot(self.start - point, d)
        t = np.dot(point - self.end, d)

        # Clamped parallel distance.
        h = np.maximum.reduce([s, t, 0])

        # Perpendicular distance component.
        c = np.cross(point - self.start, d)

        return np.hypot(h, np.linalg.norm(c))
    
@dataclass 
class Plane3D:
    """
    A plane in 3D.
    
    :param np.ndarray center: The center of the plane, with shape (3,).
    :param np.ndarray normal: The normal vector of the plane, with shape (3,).
    """
    center: np.ndarray
    normal: np.ndarray

    def same_side_of_plane(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """
        Check if two points are on the same side of a plane.
        
        :param np.ndarray p1: The first point, with shape (3,).
        :param np.ndarray p2: The second point, with shape (3,).
        :return: True if the points are on the same side of the plane, False otherwise.
        :rtype: bool
        """
        c, n = self.center, self.normal
        return np.sign(np.dot(p1 - c, n)) == np.sign(np.dot(p2 - c, n))

class CapType(Enum):
    """
    The type of a cap.
    
    :cvar NoCap: No cap.
    :cvar Flat: Flat cap.
    :cvar round: Round cap.
    """
    NoCap = auto()
    Flat = auto()
    Round = auto()

    def generate_points_on_surface(
        self, 
        center: np.ndarray, 
        radius: float,
        normal: np.ndarray, 
        res: int
    ) -> np.ndarray:
        """
        Generate points on the surface of the cap.
        
        :param np.ndarray center: The center of the cap, with shape (3,).
        :param np.ndarray normal: The normal vector of the cap, with shape (3,).
        :param int res: The resolution of the cap.
        :return: The points on the surface of the cap, with shape (N, 3).
        :rtype: np.ndarray
        """
        if self == CapType.NoCap:
            return []
        
        elif self == CapType.Flat:
            circle = Circle3D(center, radius, normal)
            return circle.generate_points_on_surface(res)
        
        elif self == CapType.Round:
            sphere = Sphere(center, radius)
            points = sphere.generate_points_on_surface(res)
            # plane = Plane3D(center, normal) 
            # TODO: Check which points are visible at all from sphere, since half
            #       of the points will be inside the cylinder anyway.
            return points
        
        else:
            raise ValueError(f"Unknown cap type: '{self}'")

@dataclass 
class Cylinder:
    start: np.ndarray
    end: np.ndarray
    radius: float
    cap_type: CapType

    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the cylinder.
        
        :param int res: The resolution of the cylinder.
        :return: The points on the surface of the cylinder, with shape (N, 3).
        :rtype: np.ndarray
        """
        res /= 4 # TODO: implement proper resolution.
        res = int(res)

        normal = self.end - self.start
        normal = normal / np.linalg.norm(normal)

        # Generate points along core of cylinder. 
        line = Line3D(self.start, self.end)
        centers = line.generate_points_along_line(res)

        # For each radius, generate points on the circumference.
        points = []

        for center in centers:
            circle = Circle3D(center, self.radius, normal)
            points_on_circle = circle.generate_points_on_circumference(res)
            points.extend(points_on_circle)

        points.extend(self.cap_type.generate_points_on_surface(self.start, self.radius, normal, res))
        points.extend(self.cap_type.generate_points_on_surface(self.end, self.radius, normal, res))

        return points
    
    def point_is_inside(self, point: np.ndarray) -> bool:
        """
        Check if a point is inside the cylinder.
        
        :param np.ndarray point: The point to check, with shape (3,).
        :return: True if the point is inside the cylinder, False otherwise.
        :rtype: bool
        """
        plane = Plane3D(point, self.end - self.start)
        if plane.same_side_of_plane(self.start, self.end):
            return False 

        else:
            line = Line3D(self.start, self.end)
            dist = line.distance_to_line(point)
            return dist <= self.radius