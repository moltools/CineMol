from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
import typing as ty

import numpy as np
    
def gram_schmidt(n: np.ndarray) -> ty.Tuple[np.ndarray, np.ndarray]:
    v = np.random.rand(3)
    v = v - np.dot(v, n) * n 
    v = v / np.linalg.norm(v)
    w = np.cross(n, v)
    return v, w

def same_side_of_plane(ref: np.ndarray, normal: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> bool:
    return np.sign(np.dot(p1 - ref, normal)) == np.sign(np.dot(p2 - ref, normal))

class Shape3D(ABC):
    @abstractmethod 
    def generate_points_on_surface(self, res: int) -> np.ndarray:
        ...

@dataclass 
class Sphere(Shape3D):
    center: np.ndarray
    radius: float

    def generate_points_on_surface(self, resolution: int) -> np.ndarray:
        phi = np.linspace(0, 2 * math.pi, resolution)
        theta = np.linspace(0, math.pi, resolution)
        phi, theta = np.meshgrid(phi, theta)
        x = self.center[0] + self.radius * np.sin(theta) * np.cos(phi)
        y = self.center[1] + self.radius * np.sin(theta) * np.sin(phi)
        z = self.center[2] + self.radius * np.cos(theta)
        m = np.array([x.flatten(), y.flatten(), z.flatten()])
        return m.T

# @dataclass
# class Circle3D(Shape3D):
#     center: np.ndarray
#     radius: float
#     normal: np.ndarray

#     def generate_points_on_circumference(self, res: int):
#         v, w = gram_schmidt(self.normal / np.linalg.norm(self.normal))
#         angles = np.linspace(0, 2 * np.pi, res)
#         return self.center + self.radius * np.outer(np.cos(angles), v) + self.radius * np.outer(np.sin(angles), w)
    
#     def generate_points_on_surface(self, resolution: int) -> np.ndarray:
#         points = []

#         radii = np.linspace(0, self.radius, resolution)
#         for radius in radii:
#             cap = Circle3D(self.center, radius, self.normal)
#             points.extend(cap.generate_points_on_circumference(resolution))

#         return points
    
#     def point_is_inside(self, _: Point3D) -> bool:
#         return False

# @dataclass 
# class Line3D(Shape3D):
#     start: Point3D
#     end: Point3D
#     radius: float

#     def generate_points_on_surface(self, _: int) -> ty.List[Point3D]:
#         return []
    
#     def point_is_inside(self, _: Point3D) -> bool:
#         return False 
    
#     def generate_points_along_line(self, resolution: int) -> ty.List[Point3D]:
#         x = np.linspace(self.start.x, self.end.x, resolution)
#         y = np.linspace(self.start.y, self.end.y, resolution)
#         z = np.linspace(self.start.z, self.end.z, resolution)
#         points = [Point3D(x[i], y[i], z[i]) for i in range(resolution)]
#         return points    

# class CapType(Enum):
#     NONE = auto()
#     FLAT = auto()
#     ROUND = auto()

# @dataclass 
# class Cylinder(Shape3D):
#     start: Point3D
#     end: Point3D
#     radius: float
#     cap_type: CapType