"""
This module contains classes for representing 3D shapes.
"""
from dataclasses import dataclass
from enum import Enum, auto
import math
import typing as ty

from cinemol.geometry import Vector3D, Point3D, sign 

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

    :param Point3D center: The center of the sphere.
    :param float radius: The radius of the sphere.
    """
    center: Point3D
    radius: float

    def generate_points_on_surface(self, res: int) -> ty.List[Point3D]:
        """
        Generate points on the surface of the sphere.
        
        :param int res: The resolution of the sphere.
        :return: The points on the surface of the sphere.
        :rtype: ty.List[Point3D]
        """
        phis, thetas = [], []
        for i in range(res):
            phis.append(2 * math.pi * i / res)
            thetas.append(math.pi * i / res)

        points = []
        for phi in phis:
            for theta in thetas:
                x = self.center.x + self.radius * math.sin(theta) * math.cos(phi)
                y = self.center.y + self.radius * math.sin(theta) * math.sin(phi)
                z = self.center.z + self.radius * math.cos(theta)

                if z > 0:
                    points.append(Point3D(x, y, z))
        
        return points
    
    def point_is_inside(self, point: Point3D) -> bool:
        """
        Check if a point is inside the sphere.
        
        :param Point3D point: The point to check.
        :return: True if the point is inside the sphere, False otherwise.
        :rtype: bool
        """
        dist = self.center.calculate_distance(point)
        return dist <= self.radius
    
@dataclass
class Circle3D:
    """
    A circle in 3D.

    :param Point3D center: The center of the circle.
    :param float radius: The radius of the circle.
    :param Vector3D normal: The normal vector of the circle.
    """
    center: Point3D
    radius: float
    normal: Vector3D    

    def generate_points_on_circumference(self, res: int) -> ty.List[Point3D]:
        """
        Generate points on the circumference of the circle.
        
        :param int res: The resolution of the circle.
        :return: The points on the circumference of the circle, with shape (N, 3).
        :rtype: ty.List[Point3D]
        """
        normal = self.normal
        normal = normal.normalize()
        v, w = gram_schmidt(normal)
        angles = []
        for i in range(res):
            angles.append(2 * math.pi * i / res)

        cos_angles = [math.cos(angle) for angle in angles]
        sin_angles = [math.sin(angle) for angle in angles]

        points = []
        for cos_angle, sin_angle in zip(cos_angles, sin_angles):
            point = Point3D(
                self.center.x + self.radius * cos_angle * v.x + self.radius * sin_angle * w.x,
                self.center.y + self.radius * cos_angle * v.y + self.radius * sin_angle * w.y,
                self.center.z + self.radius * cos_angle * v.z + self.radius * sin_angle * w.z
            )
            points.append(point)

        return points
    
    def generate_points_on_surface(self, res: int) -> ty.List[Point3D]:
        """
        Generate points on the surface of the circle.
        
        :param int res: The resolution of the circle.
        :return: The points on the surface of the circle, with shape (N, 3).
        :rtype: ty.List[Point3D]
        """
        points = []

        # For each radius, generate points on the circumference.
        radii = []
        for i in range(res):
            radii.append(self.radius * i / res)

        for radius in radii:
            cap = Circle3D(self.center, radius, self.normal)
            points.extend(cap.generate_points_on_circumference(res))

        return points

@dataclass 
class Line3D:
    """
    A line in 3D.
    
    :param Point3D start: The start point of the line.
    :param Point3D end: The end point of the line.
    """
    start: Point3D
    end: Point3D
    
    def generate_points_along_line(self, res: int) -> ty.List[Point3D]:
        """
        Generate points along the line.

        :param int res: The resolution of the line.
        :return: The points along the line.
        :rtype: ty.List[Point3D]    
        """
        points = []
        for i in range(res):
            x = self.start.x + (self.end.x - self.start.x) * i / res
            y = self.start.y + (self.end.y - self.start.y) * i / res
            z = self.start.z + (self.end.z - self.start.z) * i / res
            points.append(Point3D(x, y, z))

        return points
    
    def distance_to_line(self, point: Point3D) -> float:
        """
        Compute the distance from a point to the line.

        :param Point3D point: The point to compute the distance to.
        :return: The distance from the point to the line.
        :rtype: float
        """
        d = self.end.create_vector(self.start).normalize() # Normalized tangent
        s = self.start.create_vector(point).dot(d) # Signd parallel distance component 1
        t = point.create_vector(self.end).dot(d) # Signed parallel distance component 2
        h = max(s, t, 0.0) # Clamped parallel distance component
        c = point.create_vector(self.start).cross(d).length() # Perpendicular distance component
        hypot = math.sqrt(h * h + c * c) # Distance from point to line (Pythagoras)
        return hypot 
    
@dataclass 
class Plane3D:
    """
    A plane in 3D.
    
    :param Point3D center: The center of the plane.
    :param Vector3D normal: The normal vector of the plane.   
    """
    center: Point3D
    normal: Vector3D

    def same_side_of_plane(self, p1: Point3D, p2: Point3D) -> bool:
        """
        Check if two points are on the same side of a plane.
        
        :param Point3D p1: The first point.
        :param Point3D p2: The second point.
        :return: True if the points are on the same side of the plane, False otherwise.
        :rtype: bool
        """
        left = sign(p1.create_vector(self.center).dot(self.normal))
        right = sign(p2.create_vector(self.center).dot(self.normal))
        return left == right 

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
        center: Point3D, 
        radius: float,
        normal: Point3D, 
        res: int
    ) -> ty.List[Point3D]:
        """
        Generate points on the surface of the cap.
        
        :param Point3D center: The center of the cap.
        :param float radius: The radius of the cap.
        :param int res: The resolution of the cap.
        :return: The points on the surface of the cap.
        :rtype: ty.List[Point3D]
        """
        if self == CapType.NoCap:
            return []
        
        elif self == CapType.Flat:
            circle = Circle3D(center, radius, normal)
            return circle.generate_points_on_surface(res)
        
        elif self == CapType.Round:
            sphere = Sphere(center, radius)
            points = sphere.generate_points_on_surface(res)
            return points
        
        else:
            raise ValueError(f"Unknown cap type: '{self}'")

@dataclass 
class Cylinder:
    start: Point3D
    end: Point3D
    radius: float
    cap_type: CapType

    def generate_points_on_surface(self, res: int) -> ty.List[Point3D]:
        """
        Generate points on the surface of the cylinder.
        
        :param int res: The resolution of the cylinder.
        :return: The points on the surface of the cylinder. 
        :rtype: ty.List[Point3D]
        """
        normal = self.end.create_vector(self.start).normalize()
        centers = Line3D(self.start, self.end).generate_points_along_line(res)

        points = []
        for center in centers:
            circle = Circle3D(center, self.radius, normal)
            points_on_circle = circle.generate_points_on_circumference(res)
            points.extend(points_on_circle)

        points.extend(self.cap_type.generate_points_on_surface(self.start, self.radius, normal, res))
        points.extend(self.cap_type.generate_points_on_surface(self.end, self.radius, normal, res))

        return points
    
    def point_is_inside(self, point: Point3D) -> bool:
        """
        Check if a point is inside the cylinder.
        
        :param Point3D point: The point to check.
        :return: True if the point is inside the cylinder, False otherwise.
        :rtype: bool
        """
        plane = Plane3D(point, self.end.create_vector(self.start))
        if plane.same_side_of_plane(self.start, self.end):
            return False 

        else:
            line = Line3D(self.start, self.end)
            dist = line.distance_to_line(point)
            return dist <= self.radius