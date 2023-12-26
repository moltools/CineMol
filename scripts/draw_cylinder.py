import math 
import typing as ty
from dataclasses import dataclass
from enum import Enum 

import matplotlib.pyplot as plt
import numpy as np

class Cap(Enum):
    Round = 1
    Flat = 2
    NoCap = 3

@dataclass
class Point3D:
    x: float
    y: float
    z: float

@dataclass
class Vector3D:
    x: float
    y: float
    z: float

@dataclass
class Circle3D:
    center: Point3D
    radius: float
    normal: Vector3D

def normalize(V):
    return V / np.linalg.norm(V)

def gram_schmidt(N):
    V = np.random.rand(3)
    V = V - np.dot(V, N) * N # Create V orthogonal to N
    V = normalize(V)
    W = np.cross(N, V) # Create W orthogonal to N and V
    return V, W

def calculate_normal(p1: Point3D, p2: Point3D) -> Vector3D:
    return Vector3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)

def calculate_points_along_line(start: Point3D, end: Point3D, num_points: int) -> ty.List[Point3D]:
    x = np.linspace(start.x, end.x, num_points)
    y = np.linspace(start.y, end.y, num_points)
    z = np.linspace(start.z, end.z, num_points)
    points = [Point3D(x[i], y[i], z[i]) for i in range(num_points)]
    return points

def generate_points_on_circle(circle: Circle3D, num_points: int):
    normal_vector = np.array([circle.normal.x, circle.normal.y, circle.normal.z])
    center = np.array([circle.center.x, circle.center.y, circle.center.z])
    radius = circle.radius
    N = normalize(normal_vector)
    V, W = gram_schmidt(N) # TODO: for cylinder we calculate this now many times...
    angles = np.linspace(0, 2 * np.pi, num_points)
    points = center + radius * np.outer(np.cos(angles), V) + radius * np.outer(np.sin(angles), W)
    points = [Point3D(point[0], point[1], point[2]) for point in points]
    return points


def generate_points_on_sphere(center: Point3D, radius: float, resolution: int) -> ty.List[Point3D]:
    num_points = resolution
    sphere1_x_center = center.x
    sphere1_y_center = center.y
    sphere1_z_center = center.z
    phi = np.linspace(0, 2 * np.pi, num_points)
    theta = np.linspace(0, np.pi, num_points)
    phi, theta = np.meshgrid(phi, theta)
    x_points = sphere1_x_center + radius * np.sin(theta) * np.cos(phi)
    y_points = sphere1_y_center + radius * np.sin(theta) * np.sin(phi)
    z_points = sphere1_z_center + radius * np.cos(theta)
    points = np.array([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T
    points = [Point3D(point[0], point[1], point[2]) for point in points]
    return points

def point_side_of_plane(plane_normal, plane_point, point1, point2):
    vector1 = np.array(point1) - np.array(plane_point)
    vector2 = np.array(point2) - np.array(plane_point)
    dot_product1 = np.dot(vector1, plane_normal)
    dot_product2 = np.dot(vector2, plane_normal)
    return np.sign(dot_product1) == np.sign(dot_product2)

def calc_cap(normal: Vector3D, center: Point3D, other: Point3D, radius: float, resolution: int, cap: Cap, points: ty.List[Point3D]):
    normal = np.array([normal.x, normal.y, normal.z])
    match cap:
        case Cap.Round:
            points = generate_points_on_sphere(center, radius, resolution)
            points = np.array([[point.x, point.y, point.z] for point in points])
            plane_normal = np.array([other.x - center.x, other.y - center.y, other.z - center.z])
            plane_normal = normalize(plane_normal)
            center = np.array([center.x, center.y, center.z])
            other = np.array([other.x, other.y, other.z])
            points = [point for point in points if not point_side_of_plane(plane_normal, center, point, other)]
            # TODO: can probably speed up calc above by just calculating two opositie points on sphere and then just checking if point is on the same side as those two points
            points = [Point3D(point[0], point[1], point[2]) for point in points]
        case Cap.Flat:
            points = []
            radii = np.linspace(0, radius, resolution)
            for radius in radii:
                cap1 = Circle3D(center, radius, normal)
                points.extend(generate_points_on_circle(cap1, resolution))
        case Cap.NoCap:
            points= []
    return points 

def generate_points_on_cylinder(start: Point3D, end: Point3D, radius: float, resolution: int, start_cap: Cap, end_cap: Cap) -> ty.List[Point3D]:
    normal = calculate_normal(start, end)
    centers = calculate_points_along_line(start, end, resolution)
    points = []
    for center in centers:
        circle = Circle3D(center, radius, normal)
        points_on_circle = generate_points_on_circle(circle, resolution)
        points.extend(points_on_circle)
    
    start_cap_points = calc_cap(normal, start, end, radius, resolution, start_cap, points)
    end_cap_points = calc_cap(normal, end, start, radius, resolution, end_cap, points)
    points.extend(start_cap_points)
    points.extend(end_cap_points)
    return points
        
def plot_points(points: ty.List[Point3D], colors: ty.List[int]) -> None:
    x = [point.x for point in points]
    y = [point.y for point in points]
    z = [point.z for point in points]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
    ax.scatter(x, y, z, color=colors, s=1)
    plt.tight_layout()
    plt.show()
    plt.clf()
    
def main() -> None:
    points = generate_points_on_cylinder(Point3D(0, 0, 0), Point3D(0, 5, 5), 1, 30, Cap.Round, Cap.Round)
    sphere_center = Point3D(0, 5.5, 4)
    sphere = generate_points_on_sphere(sphere_center, 2, 30)
    colors = []
    for point in points: 
        dist = math.sqrt((point.x - sphere_center.x) ** 2 + (point.y - sphere_center.y) ** 2 + (point.z - sphere_center.z) ** 2)
        if dist < 2: colors.append("red")
        else: colors.append("blue")
    points += sphere
    for point in sphere: colors.append("green")
    plot_points(points, colors)
    exit(0)

if __name__ == "__main__":
    main()
