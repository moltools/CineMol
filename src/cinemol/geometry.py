from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto    
import math
import random 
import typing as ty

import numpy as np # TODO: remove numpy dependency

# ==============================================================================
# Graph geometry
# ==============================================================================

@dataclass
class Point2D:
    x: float 
    y: float 

@dataclass
class Vector3D:
    x: float
    y: float
    z: float

    @classmethod 
    def random(cls) -> "Vector3D":
        return Vector3D(random.random(), random.random(), random.random())

    def normalize(self) -> "Vector3D":
        length = (self.x ** 2 + self.y ** 2 + self.z ** 2) ** 0.5
        return Vector3D(self.x / length, self.y / length, self.z / length)
    
    def dot(self, other: "Vector3D") -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z
    
    def cross(self, other: "Vector3D") -> "Vector3D":
        return Vector3D(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z, 
            self.x * other.y - self.y * other.x
        )
    
    def gram_schmidt(self) -> ty.Tuple["Vector3D", "Vector3D"]:
        v1 = Vector3D.random()
        v2 = v1.dot(self)
        v3 = Vector3D(v1.x - v2 * self.x, v1.y - v2 * self.y, v1.z - v2 * self.z).normalize()
        v4 = self.cross(v3)
        return v3, v4

@dataclass
class Point3D:
    x: float
    y: float
    z: float

    def create_vector(self, other: "Point3D") -> Vector3D:
        return Vector3D(other.x - self.x, other.y - self.y, other.z - self.z)

def calculate_distance(a: Point3D, b: Point3D) -> float:
    return ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z)) ** 0.5

# ==============================================================================
# 3D Shapes
# ==============================================================================

def same_side_of_plane(plane_point: Point3D, plane_normal: Vector3D, p1: Point3D, p2: Point3D) -> bool:
    plane_point = np.array([plane_point.x, plane_point.y, plane_point.z])
    plane_normal = np.array([plane_normal.x, plane_normal.y, plane_normal.z])
    p1 = np.array([p1.x, p1.y, p1.z])
    p2 = np.array([p2.x, p2.y, p2.z])

    vector1 = np.array(p1) - np.array(plane_point)
    vector2 = np.array(p2) - np.array(plane_point)
    dot_product1 = np.dot(vector1, plane_normal)
    dot_product2 = np.dot(vector2, plane_normal)
    return np.sign(dot_product1) == np.sign(dot_product2)

class Shape3D(ABC):
    @abstractmethod 
    def generate_points_on_surface(self, reslution: int) -> ty.List[Point3D]:
        ...

    @abstractmethod
    def point_is_inside(self, point: Point3D) -> str:
        ...

@dataclass 
class Sphere(Shape3D):
    center: Point3D
    radius: float

    def generate_points_on_surface(self, resolution: int) -> ty.List[Point3D]:
        phi = np.linspace(0, 2 * math.pi, resolution)
        theta = np.linspace(0, math.pi, resolution)
        phi, theta = np.meshgrid(phi, theta)
        x_points = self.center.x + self.radius * math.sin(theta) * math.cos(phi)
        y_points = self.center.y + self.radius * math.sin(theta) * math.sin(phi)
        z_points = self.center.z + self.radius * math.cos(theta)
        points = np.array([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T
        points = [Point3D(point[0], point[1], point[2]) for point in points]
        return points
    
    def point_is_inside(self, point: Point3D) -> bool:
        return calculate_distance(self.center, point) <= self.radius

@dataclass
class Circle3D(Shape3D):
    center: Point3D
    radius: float
    normal: Vector3D

    def generate_points_on_circumference(self, resolution: int):
        v1 = Point3D(0, 0, 0).create_vector(self.normal).normalize()
        v2, v3 = v1.gram_schmidt()
        angles = np.linspace(0, 2 * np.pi, resolution)
        
        v2 = np.array([v2.x, v2.y, v2.z])
        v3 = np.array([v3.x, v3.y, v3.z])
        center = np.array([self.center.x, self.center.y, self.center.z])
        points = center + self.radius * np.outer(np.cos(angles), v2) + self.radius * np.outer(np.sin(angles), v3)
        
        points = [Point3D(point[0], point[1], point[2]) for point in points]
        return points
    
    def generate_points_on_surface(self, resolution: int) -> ty.List[Point3D]:
        points = []

        radii = np.linspace(0, self.radius, resolution)
        for radius in radii:
            cap = Circle3D(self.center, radius, self.normal)
            points.extend(cap.generate_points_on_circumference(resolution))

        return points
    
    def point_is_inside(self, _: Point3D) -> bool:
        return False

@dataclass 
class Line3D(Shape3D):
    start: Point3D
    end: Point3D
    radius: float

    def generate_points_on_surface(self, _: int) -> ty.List[Point3D]:
        return []
    
    def point_is_inside(self, _: Point3D) -> bool:
        return False 
    
    def generate_points_along_line(self, resolution: int) -> ty.List[Point3D]:
        x = np.linspace(self.start.x, self.end.x, resolution)
        y = np.linspace(self.start.y, self.end.y, resolution)
        z = np.linspace(self.start.z, self.end.z, resolution)
        points = [Point3D(x[i], y[i], z[i]) for i in range(resolution)]
        return points    

class CapType(Enum):
    NONE = auto()
    FLAT = auto()
    ROUND = auto()

@dataclass 
class Cylinder(Shape3D):
    start: Point3D
    end: Point3D
    radius: float
    cap_type: CapType

# ==============================================================================
# 2D Shapes
# ==============================================================================

class Shape2D(ABC):
    @abstractmethod
    def to_svg(self) -> str:
        ...

@dataclass 
class Circle2D(Shape2D):
    reference: str
    center: Point2D
    radius: float   

    def to_svg(self) -> str:
        cx, cy, r = self.center.x, self.center.y, self.radius
        return f'<circle class="{self.reference}" cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}"/>'

@dataclass
class Line2D(Shape2D):
    reference: str
    start: Point2D
    end: Point2D

    def to_svg(self) -> str:
        x1, y1, x2, y2 = self.start.x, self.start.y, self.end.x, self.end.y
        return f'<line class="{self.reference}" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'

@dataclass
class Polygon2D(Shape2D):
    reference: str 
    points: ty.List[Point2D] # Ordered list of points

    def to_svg(self) -> str:
        points = " ".join(f"{p.x:.3f},{p.y:.3f}" for p in self.points)
        return f'<polygon class="{self.reference}" points="{points}"/>'