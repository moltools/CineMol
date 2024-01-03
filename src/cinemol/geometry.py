"""
This module contains classes for representing 3D shapes.
"""
from dataclasses import dataclass
from enum import Enum, auto
import math
import random
import typing as ty

# ==============================================================================
# Basic geometry classes
# ==============================================================================

class Vector3D:
    """
    Represents a vector in 3D space.
    """
    def __init__(self, x: float, y: float, z: float) -> None:
        """
        :param float x: x-coordinate of the vector.
        :param float y: y-coordinate of the vector.
        :param float z: z-coordinate of the vector.
        """
        self.x = x 
        self.y = y
        self.z = z

    @classmethod
    def create_random(cls) -> "Vector3D":
        """
        Creates a random vector.
        
        :return: A random vector.
        :rtype: Vector3D
        """
        return Vector3D(random.random(), random.random(), random.random())
    
    def length(self) -> float:
        """
        Calculates the length of this vector.
        
        :return: The length of this vector.
        :rtype: float
        """
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
    
    def normalize(self) -> "Vector3D":
        """
        Normalizes this vector.
        
        :return: A new vector with the same direction as this vector but with unit length.
        :rtype: Vector3D
        """
        return self.multiply(1 / self.length())
    
    def dot(self, other: "Vector3D") -> float:
        """
        Calculates the dot product of this vector with another vector.
        
        :param Vector3D other: The vector to dot with this vector.
        :return: The dot product of this vector with another vector.
        :rtype: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z
    
    def cross(self, other: "Vector3D") -> "Vector3D":
        """
        Calculates the cross product of this vector with another vector.
        
        :param Vector3D other: The vector to cross with this vector.
        :return: The cross product of this vector with another vector.
        :rtype: Vector3D
        """
        return Vector3D(
            self.y * other.z - self.z * other.y, 
            self.z * other.x - self.x * other.z, 
            self.x * other.y - self.y * other.x
        )
    
    def subtract(self, other: "Vector3D") -> "Vector3D":
        """
        Subtracts another vector from this vector.
        
        :param Vector3D other: The vector to subtract from this vector.
        :return: A new vector with the coordinates of the difference.
        :rtype: Vector3D
        """
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def multiply(self, scalar: float) -> "Vector3D":
        """
        Multiplies this vector by a scalar.
        
        :param float scalar: The scalar to multiply this vector by.
        :return: A new vector with the coordinates of the product.
        :rtype: Vector3D
        """
        return Vector3D(self.x * scalar, self.y * scalar, self.z * scalar)

class Point2D:
    """
    Represents a point in 2D space.
    """
    def __init__(self, x: float, y: float) -> None:
        """
        :param float x: x-coordinate of the point.
        :param float y: y-coordinate of the point.
        """
        self.x = x 
        self.y = y

    def subtract_point(self, other: "Point2D") -> "Point2D":
        """
        Subtracts the coordinates of another point from this point.
        
        :param Point2D other: The point to subtract from this point.
        :return: A new point with the coordinates of the difference.
        :rtype: Point2D
        """
        return Point2D(self.x - other.x, self.y - other.y)
    
    def cross(self, other: "Point2D") -> float:
        """
        Calculates the cross product of this point with another point.
        
        :param Point2D other: The point to cross with this point.
        :return: The cross product of this point with another point.
        :rtype: float
        """
        return self.x * other.y - self.y * other.x
    
class Point3D:
    """
    Represents a point in 3D space.
    """
    def __init__(self, x: float, y: float, z: float) -> None:
        """
        :param float x: x-coordinate of the point.
        :param float y: y-coordinate of the point.
        :param float z: z-coordinate of the point.
        """
        self.x = x 
        self.y = y
        self.z = z
    
    def create_vector(self, other: "Point3D") -> Vector3D:
        """
        Creates a vector from this point to another point.
        
        :param Point3D other: The other point.
        :return: A vector from this point to another point.
        :rtype: Vector3D
        """
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def calculate_distance(self, other: "Point3D") -> float:
        """
        Calculates the distance between this point and another point.

        :param Point3D other: The other point.
        :return: The distance between this point and another point.
        :rtype: float
        """
        return self.create_vector(other).length()

    def midpoint(self, other: "Point3D") -> "Point3D":
        """
        Calculates the midpoint between this point and another point.

        :param Point3D other: The other point.
        :return: The midpoint between this point and another point.
        :rtype: Point3D
        """
        return Point3D((self.x + other.x) / 2, (self.y + other.y) / 2, (self.z + other.z) / 2)
    
    def translate(self, vector: "Vector3D") -> "Point3D":
        """
        Adds a vector to this point.
        
        :param Vector3D vector: The vector to add to this point.
        :return: A new point with the coordinates of the sum.
        :rtype: Point3D
        """
        return Point3D(self.x + vector.x, self.y + vector.y, self.z + vector.z)
    
    def rotate(
        self, 
        x: float = 0.0, 
        y: float = 0.0, 
        z: float = 0.0
    ) -> "Point3D":
        """
        Rotates this point around the origin.
        
        :param float x: The clockwise rotation around the x-axis.
        :param float y: The clockwise rotation around the y-axis.
        :param float z: The clockwise rotation around the z-axis.
        :return: A new point with the coordinates of the rotated point.
        :rtype: Point3D

        Note x, y, z and in radians.
        """
        # Rotate around x-axis.
        y1 = self.y * math.cos(x) - self.z * math.sin(x)
        z1 = self.y * math.sin(x) + self.z * math.cos(x)

        # Rotate around y-axis.
        x2 = self.x * math.cos(y) + z1 * math.sin(y)
        z2 = -self.x * math.sin(y) + z1 * math.cos(y)

        # Rotate around z-axis.
        x3 = x2 * math.cos(z) - y1 * math.sin(z)
        y3 = x2 * math.sin(z) + y1 * math.cos(z)

        return Point3D(x3, y3, z2)

# ==============================================================================
# Helper functions
# ==============================================================================

def sign(x: float) -> int:
    """
    Returns the sign of a number.
    
    :param float x: The number.
    :return: The sign of the number.
    :rtype: int
    """
    if x < 0:
        return -1
    elif x > 0:
        return 1
    else:
        return 0

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
# Shape definitions
# ==============================================================================

@dataclass 
class Line3D:
    """
    A line in 3D.
    
    :param Point3D start: The start point of the line.
    :param Point3D end: The end point of the line.
    """
    start: Point3D
    end: Point3D

@dataclass 
class Plane3D:
    """
    A plane in 3D.
    
    :param Point3D center: The center of the plane.
    :param Vector3D normal: The normal vector of the plane.   
    """
    center: Point3D
    normal: Vector3D

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

@dataclass 
class Sphere:
    """
    A sphere.

    :param Point3D center: The center of the sphere.
    :param float radius: The radius of the sphere.
    """
    center: Point3D
    radius: float

class CylinderCapType(Enum):
    """
    The type of a cap.
    
    :cvar NoCap: No cap.
    :cvar Flat: Flat cap.
    :cvar round: Round cap.
    """
    NoCap = auto()
    Flat = auto()
    Round = auto()

@dataclass 
class Cylinder:
    start: Point3D
    end: Point3D
    radius: float
    cap_type: CylinderCapType

# ==============================================================================
# Check if points are on the same side of a plane
# ============================================================================== 

def same_side_of_plane(plane: Plane3D, p1: Point3D, p2: Point3D) -> bool:
    """
    Check if two points are on the same side of a plane.
    
    :param Plane3D plane: The plane.
    :param Point3D p1: The first point.
    :param Point3D p2: The second point.
    :return: True if the points are on the same side of the plane, False otherwise.
    :rtype: bool
    """
    left = sign(p1.create_vector(plane.center).dot(plane.normal))
    right = sign(p2.create_vector(plane.center).dot(plane.normal))
    return left == right 

# ==============================================================================
# Compute distance from point to line
# ============================================================================== 

def distance_to_line(line: Line3D, point: Point3D) -> float:
    """
    Compute the distance from a point to the line.

    :param Line3D line: The line to compute the distance to.
    :param Point3D point: The point to compute the distance to.
    :return: The distance from the point to the line.
    :rtype: float
    """
    d = line.end.create_vector(line.start).normalize()
    s = line.start.create_vector(point).dot(d) 
    t = point.create_vector(line.end).dot(d) 
    h = max(s, t, 0.0) 
    c = point.create_vector(line.start).cross(d).length() 
    return math.sqrt(h * h + c * c) # Pythagorean theorem.

# ==============================================================================
# Get perpendicular lines
# ============================================================================== 

def get_perpendicular_lines(line: Line3D, width: float, num_lines: int) -> ty.List[Line3D]:
    """
    Split current line into multiple perpendicular lines.

    :param Line3D line: The line to split.
    :param float width: The spacing between the lines.
    :param int num_lines: The number of lines to split the line into.
    :return: The perpendicular lines.
    :rtype: ty.List[Line3D]
    """
    if num_lines == 1:
        return [line]

    # Ignore z-axis and get a vector perpendicular to the line to translate start and end points on.
    v = Vector3D(line.end.y - line.start.y, line.start.x - line.end.x, 0.0).normalize()

    # Get new start and end points, but centroid of starts and ends should always be original start and end.
    start = line.start.translate(v.multiply(-width * (num_lines - 1) / 2))
    end = line.end.translate(v.multiply(-width * (num_lines - 1) / 2))

    # Get new lines.
    lines = []
    for _ in range(num_lines):
        lines.append(Line3D(start, end))
        start = start.translate(v.multiply(width))
        end = end.translate(v.multiply(width))
    
    return lines
    

# ==============================================================================
# Generate points based on shape
# ==============================================================================
    
def get_points_on_line_3d(line: Line3D, num_points: int) -> ty.List[Point3D]:
    """
    Generate points along a line.

    :param Line3D line: The line.
    :param int num_points: The number of points to generate.
    :return: The points along the line.
    :rtype: ty.List[Point3D]
    """
    s_cx, s_cy, s_cz = line.start.x, line.start.y, line.start.z
    e_cx, e_cy, e_cz = line.end.x, line.end.y, line.end.z

    points = []
    for i in range(num_points + 1):
        i /= num_points 
        point = Point3D(
            s_cx + (e_cx - s_cx) * i,
            s_cy + (e_cy - s_cy) * i,
            s_cz + (e_cz - s_cz) * i
        )
        points.append(point)

    return points

def get_points_on_circumference_circle_3d(circle: Circle3D, num_points: int) -> ty.List[Point3D]:
    """
    Generate points on the circumference of the circle.
    
    :param int res: The resolution of the circle.
    :return: The points on the circumference of the circle, with shape (N, 3).
    :rtype: ty.List[Point3D]
    """
    angles = [2 * math.pi * i / num_points for i in range(num_points)]
    normal = circle.normal.normalize()
    v, w = gram_schmidt(normal)

    cos_angles = [math.cos(angle) for angle in angles]
    sin_angles = [math.sin(angle) for angle in angles]

    cx, cy, cz = circle.center.x, circle.center.y, circle.center.z
    r = circle.radius

    points = []
    for cos_a, sin_a in zip(cos_angles, sin_angles):
        point = Point3D(
            cx + r * cos_a * v.x + r * sin_a * w.x,
            cy + r * cos_a * v.y + r * sin_a * w.y,
            cz + r * cos_a * v.z + r * sin_a * w.z
        )
        points.append(point)

    return points

def get_points_on_surface_circle_3d(
    circle: Circle3D, 
    num_radii: int, 
    num_points: int
) -> ty.List[Point3D]:
    """
    Generate points on the surface of the circle.

    :param Circle3D circle: The circle.
    :param int num_radii: The number of radii to generate points for.
    :param int num_points: The number of points to generate per radius.
    :return: The points on the surface of the circle.
    :rtype: ty.List[Point3D]
    """
    radii = [circle.radius * i / num_radii for i in range(num_radii)]

    points = []
    for radius in radii:
        temp_circle = Circle3D(circle.center, radius, circle.normal)
        points.extend(get_points_on_circumference_circle_3d(temp_circle, num_points))

    return points

def get_points_on_surface_sphere(
    sphere: Sphere, 
    num_phi: int, 
    num_theta: int,
    filter_for_pov: bool = True
) -> ty.List[Point3D]:
    """
    Generate points on the surface of a sphere.
    
    :param Sphere sphere: The sphere.
    :param int num_phi: The resolution of the sphere in the phi direction.
    :param int num_theta: The resolution of the sphere in the theta direction.
    :param ty.Optional[Vector3D] filter_for_pov: If True, only return points 
        that are on the surface of the sphere we can see from POV positive z-axis towards origin.
    :return: The points on the surface of the sphere.
    :rtype: ty.List[Point3D]
    """
    phis = [2 * math.pi * i / num_phi for i in range(num_phi + 1)]
    thetas = [math.pi * i / num_theta for i in range(num_theta + 1)]

    cx, cy, cz = sphere.center.x, sphere.center.y, sphere.center.z
    r = sphere.radius
    
    points = []
    for theta in thetas:
        for phi in phis:
            x = cx + r * math.sin(theta) * math.cos(phi)
            y = cy + r * math.sin(theta) * math.sin(phi)
            z = cz + r * math.cos(theta)

            # Only add points that are on the surface of the sphere we can see.
            if not filter_for_pov:
                points.append(Point3D(x, y, z))
                continue

            # Check if point is on the surface of the sphere we can see.
            if z >= sphere.center.z:
                points.append(Point3D(x, y, z))
    
    return points

def get_points_on_surface_cap(
    cap_type: CylinderCapType,
    center_cap: Point3D, 
    radius_cap: float,
    normal_cap: Point3D,
    center_cylinder: Point3D, 
    resolution: int,
    filter_for_pov: bool = True
) -> ty.List[Point3D]:
    """
    Generate points on the surface of the cap.
    
    :param CylinderCapType cap_type: The type of the cap.
    :param Point3D center_cap: The center of the cap.
    :param float radius_cap: The radius of the cap.
    :param Vector3D normal_cap: The normal vector of the cap.
    :param Point3D center_cylinder: The center of the cylinder.
    :param int resolution: The resolution of the cap.
    :param bool filter_for_pov: If True, only return points that are on the visible surface of the cap
        (i.e. the surface of the cap we can see from POV positive z-axis towards origin).
    :return: The points on the surface of the cap.
    :rtype: ty.List[Point3D]
    """
    if cap_type == CylinderCapType.NoCap:
        return []
    
    elif cap_type == CylinderCapType.Flat:
        circle = Circle3D(center_cap, radius_cap, normal_cap)
        return get_points_on_circumference_circle_3d(circle, resolution)
    
    elif cap_type == CylinderCapType.Round:
        sphere = Sphere(center_cap, radius_cap)
        plane = Plane3D(center_cap, normal_cap)
        points = get_points_on_surface_sphere(sphere, resolution, resolution, filter_for_pov=filter_for_pov)
        points = [point for point in points if not same_side_of_plane(plane, center_cylinder, point)]
        return points
    
    else:
        raise ValueError(f"Unknown cap type: '{cap_type}'")

def get_points_on_surface_cylinder(cylinder: Cylinder, resolution: int,) -> ty.List[Point3D]:
    """
    Generate points on the surface of the cylinder.
    
    :param Cylinder cylinder: The cylinder.
    :param int resolution: The resolution of the cylinder.
    :return: The points on the surface of the cylinder. 
    :rtype: ty.List[Point3D]
    """
    normal = cylinder.end.create_vector(cylinder.start).normalize()
    centers = get_points_on_line_3d(Line3D(cylinder.start, cylinder.end), resolution)

    points = []
    for center in centers:
        circle = Circle3D(center, cylinder.radius, normal)
        points.extend(get_points_on_circumference_circle_3d(circle, resolution))

    # Get points on the caps.
    cap_type = cylinder.cap_type
    cap_points = get_points_on_surface_cap(cap_type, cylinder.start, cylinder.radius, normal, cylinder.end, resolution, False)
    points.extend(cap_points)
    cap_points = get_points_on_surface_cap(cap_type, cylinder.end, cylinder.radius, normal, cylinder.start, resolution, False)
    points.extend(cap_points)

    return points

# ==============================================================================
# Check if points are inside a shape
# ==============================================================================

def point_is_inside_sphere(sphere: Sphere, point: Point3D) -> bool:
    """
    Check if a point is inside the sphere.
    
    :param Sphere sphere: The sphere to check.
    :param Point3D point: The point to check.
    :return: True if the point is inside the sphere, False otherwise.
    :rtype: bool
    """
    return sphere.center.calculate_distance(point) <= sphere.radius

def point_is_inside_cylinder(cylinder: Cylinder, point: Point3D) -> bool:
    """
    Check if a point is inside the cylinder.
    
    :param Cylinder cylinder: The cylinder to check.
    :param Point3D point: The point to check.
    :return: True if the point is inside the cylinder, False otherwise.
    :rtype: bool
    """
    line = Line3D(cylinder.start, cylinder.end)
    dist = distance_to_line(line, point)
    cap_type = cylinder.cap_type

    if cap_type == CylinderCapType.Round:
        return dist <= cylinder.radius
    
    elif cap_type == CylinderCapType.Flat or cap_type == CylinderCapType.NoCap:
        normal = cylinder.end.create_vector(cylinder.start).normalize()
        plane1 = Plane3D(cylinder.start, normal)
        plane2 = Plane3D(cylinder.end, normal)
        is_between_planes = (
            same_side_of_plane(plane1, point, cylinder.end) and 
            same_side_of_plane(plane2, point, cylinder.start)
        )
        return dist <= cylinder.radius and is_between_planes
    
    else:
        raise ValueError(f"Unknown cap type: '{cap_type}'")
    
# ==============================================================================
# Check if shapes intersect
# ==============================================================================
    
def sphere_intersects_with_sphere(sphere1: Sphere, sphere2: Sphere) -> bool:
    """
    Check if two spheres intersect.

    :param Sphere sphere1: The first sphere.
    :param Sphere sphere2: The second sphere.
    :return: True if the spheres intersect, False otherwise.
    :rtype: bool
    """
    c1, r1 = sphere1.center, sphere1.radius
    c2, r2 = sphere2.center, sphere2.radius
    return c1.calculate_distance(c2) <= r1 + r2

def sphere_intersects_with_cylinder(sphere: Sphere, cylinder: Cylinder) -> bool:
    """
    Check if a sphere intersects with a cylinder.
    
    :param Sphere sphere: The sphere.
    :param Cylinder cylinder: The cylinder.
    :return: True if the sphere intersects with the cylinder, False otherwise.
    :rtype: bool
    """
    d = sphere.radius + cylinder.radius
    return (
        sphere.center.calculate_distance(cylinder.start) <= d or 
        sphere.center.calculate_distance(cylinder.end) <= d
    )

def cylinder_intersects_with_cylinder(cylinder1: Cylinder, cylinder2: Cylinder) -> bool:
    """
    Check if two cylinders intersect.

    :param Cylinder cylinder1: The first cylinder.
    :param Cylinder cylinder2: The second cylinder.
    :return: True if the cylinders intersect, False otherwise.
    :rtype: bool
    """
    # TODO: currently a hack to check if start/end is same, but not really checking 
    #       if they intersect at any other position along the cylinder's direction.
    
    d = cylinder1.radius + cylinder2.radius
    return (
        cylinder1.start.calculate_distance(cylinder2.start) <= d or 
        cylinder1.start.calculate_distance(cylinder2.end) <= d or 
        cylinder1.end.calculate_distance(cylinder2.start) <= d or
        cylinder1.end.calculate_distance(cylinder2.end) <= d
    )
