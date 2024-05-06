# -*- coding: utf-8 -*-

"""This module contains classes for representing 3D shapes."""

import math
import typing as ty
from dataclasses import dataclass
from enum import Enum, auto
from random import SystemRandom

# ==============================================================================
# Basic geometry classes
# ==============================================================================


class Vector3D:
    """Represents a vector in 3D space."""

    def __init__(self, x: float, y: float, z: float) -> None:
        """Initialize a new vector.

        :param x: The x-coordinate of the vector.
        :type x: float
        :param y: The y-coordinate of the vector.
        :type y: float
        :param z: The z-coordinate of the vector.
        :type z: float
        """
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def create_random(cls) -> "Vector3D":
        """Create a random vector.

        :return: A random vector.
        :rtype: Vector3D
        """
        cryptogen = SystemRandom()
        return Vector3D(cryptogen.random(), cryptogen.random(), cryptogen.random())

    def length(self) -> float:
        """Calculate the length of this vector.

        :return: The length of this vector.
        :rtype: float
        """
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self) -> "Vector3D":
        """Normalize this vector.

        :return: A new vector with the same direction as this vector but with unit length.
        :rtype: Vector3D
        """
        if self.length() == 0:  # Prevent division by zero.
            return Vector3D(0, 0, 0)
        return self.multiply(1 / self.length())

    def dot(self, other: "Vector3D") -> float:
        """Calculate the dot product of this vector with another vector.

        :param other: The vector to dot with this vector.
        :type other: Vector3D
        :return: The dot product of this vector with another vector.
        :rtype: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: "Vector3D") -> "Vector3D":
        """Calculate the cross product of this vector with another vector.

        :param other: The vector to cross with this vector.
        :type other: Vector3D
        :return: The cross product of this vector with another vector.
        :rtype: Vector3D
        """
        return Vector3D(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def subtract(self, other: "Vector3D") -> "Vector3D":
        """Subtracts another vector from this vector.

        :param other: The vector to subtract from this vector.
        :type other: Vector3D
        :return: A new vector with the coordinates of the difference.
        :rtype: Vector3D
        """
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)

    def multiply(self, scalar: float) -> "Vector3D":
        """Multiplies this vector by a scalar.

        :param scalar: The scalar to multiply this vector by.
        :type scalar: float
        :return: A new vector with the coordinates of the product.
        :rtype: Vector3D
        """
        return Vector3D(self.x * scalar, self.y * scalar, self.z * scalar)


class Point2D:
    """Represents a point in 2D space."""

    def __init__(self, x: float, y: float) -> None:
        """Initialize a new point.

        :param x: The x-coordinate of the point.
        :type x: float
        :param y: The y-coordinate of the point.
        :type y: float
        """
        self.x = x
        self.y = y

    def subtract_point(self, other: "Point2D") -> "Point2D":
        """Subtracts the coordinates of another point from this point.

        :param other: The point to subtract from this point.
        :type other: Point2D
        :return: A new point with the coordinates of the difference.
        :rtype: Point2D
        """
        return Point2D(self.x - other.x, self.y - other.y)

    def cross(self, other: "Point2D") -> float:
        """Calculate the cross product of this point with another point.

        :param other: The point to cross with this point.
        :type other: Point2D
        :return: The cross product of this point with another point.
        :rtype: float
        """
        return self.x * other.y - self.y * other.x


class Point3D:
    """Represents a point in 3D space."""

    def __init__(self, x: float, y: float, z: float) -> None:
        """Initialize a new point.

        :param x: The x-coordinate of the point.
        :type x: float
        :param y: The y-coordinate of the point.
        :type y: float
        :param z: The z-coordinate of the point.
        :type z: float
        """
        self.x = x
        self.y = y
        self.z = z

    def create_vector(self, other: "Point3D") -> Vector3D:
        """Create a vector from this point to another point.

        :param other: The other point.
        :type other: Point3D
        :return: A vector from this point to another point.
        :rtype: Vector3D
        """
        return Vector3D(other.x - self.x, other.y - self.y, other.z - self.z)

    def calculate_distance(self, other: "Point3D") -> float:
        """Calculate the distance between this point and another point.

        :param other: The other point.
        :type other: Point3D
        :return: The distance between this point and another point.
        :rtype: float
        """
        return self.create_vector(other).length()

    def midpoint(self, other: "Point3D") -> "Point3D":
        """Calculate the midpoint between this point and another point.

        :param other: The other point.
        :type other: Point3D
        :return: The midpoint between this point and another point.
        :rtype: Point3D
        """
        return Point3D((self.x + other.x) / 2, (self.y + other.y) / 2, (self.z + other.z) / 2)

    def translate(self, vector: "Vector3D") -> "Point3D":
        """Add a vector to this point.

        :param vector: The vector to add to this point.
        :type vector: Vector3D
        :return: A new point with the coordinates of the sum.
        :rtype: Point3D
        """
        return Point3D(self.x + vector.x, self.y + vector.y, self.z + vector.z)

    def rotate(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> "Point3D":
        """Rotate this point around the origin.

        :param x: The clockwise rotation around the x-axis.
        :type x: float
        :param y: The clockwise rotation around the y-axis.
        :type y: float
        :param z: The clockwise rotation around the z-axis.
        :type z: float
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
    """Return the sign of a number.

    :param x: The number.
    :type x: float
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
    """Generate two orthogonal vectors for a given vector using the Gram-Schmidt process.

    :param n: The vector to generate orthogonal vectors for.
    :type n: Vector3D
    :return: Two orthogonal vectors.
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
    """A line in 3D.

    :param start: The start point of the line.
    :type start: Point3D
    :param end: The end point of the line.
    :type end: Point3D
    """

    start: Point3D
    end: Point3D


@dataclass
class Plane3D:
    """A plane in 3D.

    :param center: The center of the plane.
    :type center: Point3D
    :param normal: The normal vector of the plane.
    :type normal: Vector3D
    """

    center: Point3D
    normal: Vector3D


@dataclass
class Circle3D:
    """A circle in 3D.

    :param center: The center of the circle.
    :type center: Point3D
    :param radius: The radius of the circle.
    :type radius: float
    :param normal: The normal vector of the circle.
    :type normal: Vector3D
    """

    center: Point3D
    radius: float
    normal: Vector3D


@dataclass
class Sphere:
    """A sphere.

    :param center: The center of the sphere.
    :type center: Point3D
    :param radius: The radius of the sphere.
    :type radius: float
    """

    center: Point3D
    radius: float


class CylinderCapType(Enum):
    """
    The type of a cap.

    :cvar NO_CAP: No cap.
    :cvar FLAT: Flat cap.
    :cvar ROUND: Round cap.
    """

    NO_CAP = auto()
    FLAT = auto()
    ROUND = auto()


@dataclass
class Cylinder:
    """A cylinder.

    :param start: The start point of the cylinder.
    :type start: Point3D
    :param end: The end point of the cylinder.
    :type end: Point3D
    :param radius: The radius of the cylinder.
    :type radius: float
    :param cap_type: The type of the cap.
    :type cap_type: CylinderCapType
    """

    start: Point3D
    end: Point3D
    radius: float
    cap_type: CylinderCapType


# ==============================================================================
# Check if points are on the same side of a plane
# ==============================================================================


def same_side_of_plane(plane: Plane3D, p1: Point3D, p2: Point3D) -> bool:
    """Check if two points are on the same side of a plane.

    :param plane: The plane.
    :type plane: Plane3D
    :param p1: The first point.
    :type p1: Point3D
    :param p2: The second point.
    :type p2: Point3D
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
    """Compute the distance from a point to the line.

    :param line: The line to compute the distance to.
    :type line: Line3D
    :param point: The point to compute the distance to.
    :type point: Point3D
    :return: The distance from the point to the line.
    :rtype: float
    """
    d = line.end.create_vector(line.start).normalize()
    s = line.start.create_vector(point).dot(d)
    t = point.create_vector(line.end).dot(d)
    h = max(s, t, 0.0)
    c = point.create_vector(line.start).cross(d).length()
    return math.sqrt(h * h + c * c)  # Pythagorean theorem.


# ==============================================================================
# Get perpendicular lines
# ==============================================================================


def get_perpendicular_lines(line: Line3D, width: float, num_lines: int) -> ty.List[Line3D]:
    """Split current line into multiple perpendicular lines.

    :param line: The line to split.
    :type line: Line3D
    :param width: The spacing between the lines.
    :type width: float
    :param num_lines: The number of lines to split the line into.
    :type num_lines: int
    :return: The perpendicular lines.
    :rtype: ty.List[Line3D]
    :raises ValueError: If the number of lines is less than 1.
    :raises ValueError: If the width is less than or equal to 0.
    """
    if num_lines < 1:
        raise ValueError("Number of lines must be greater than 0.")

    if width <= 0:
        raise ValueError("Width must be greater than 0.")

    if num_lines == 1:
        return [line]

    # Ignore z-axis and get a vector perpendicular to the line to translate start
    # and end points on.
    v = Vector3D(line.end.y - line.start.y, line.start.x - line.end.x, 0.0).normalize()

    # Get new start and end points, but centroid of starts and ends should always
    # be original start and end.
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
    """Generate `num_points` + 1 points along a line.

    :param line: The line.
    :type line: Line3D
    :param num_points: The number of points to generate.
    :type num_points: int
    :return: The points along the line.
    :rtype: ty.List[Point3D]
    """
    s_cx, s_cy, s_cz = line.start.x, line.start.y, line.start.z
    e_cx, e_cy, e_cz = line.end.x, line.end.y, line.end.z

    points = []
    for i in range(num_points + 1):
        interpolation_factor: float = i / num_points
        point = Point3D(
            s_cx + (e_cx - s_cx) * interpolation_factor,
            s_cy + (e_cy - s_cy) * interpolation_factor,
            s_cz + (e_cz - s_cz) * interpolation_factor,
        )
        points.append(point)

    return points


def get_points_on_circumference_circle_3d(circle: Circle3D, num_points: int) -> ty.List[Point3D]:
    """Generate points on the circumference of the circle.

    :param circle: The circle.
    :type circle: Circle3D
    :param num_points: The number of points to generate.
    :type num_points: int
    :return: The points on the circumference of the circle.
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
            cz + r * cos_a * v.z + r * sin_a * w.z,
        )
        points.append(point)

    return points


def get_points_on_surface_circle_3d(
    circle: Circle3D, num_radii: int, num_points: int
) -> ty.List[Point3D]:
    """Generate points on the surface of the circle.

    :param circle: The circle.
    :type circle: Circle3D
    :param num_radii: The number of radii to generate points for.
    :type num_radii: int
    :param num_points: The number of points to generate per radius.
    :type num_points: int
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
    sphere: Sphere, num_phi: int, num_theta: int, filter_for_pov: bool = True
) -> ty.List[Point3D]:
    """Generate points on the surface of a sphere.

    :param sphere: The sphere.
    :type sphere: Sphere
    :param num_phi: The resolution of the sphere in the phi direction.
    :type num_phi: int
    :param num_theta: The resolution of the sphere in the theta direction.
    :type num_theta: int
    :param filter_for_pov: If True, only return points that are on the surface of
        the sphere we can see from POV positive z-axis towards origin.
    :type filter_for_pov: bool
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
    normal_cap: Vector3D,
    center_cylinder: Point3D,
    resolution: int,
    filter_for_pov: bool = True,
) -> ty.List[Point3D]:
    """Generate points on the surface of the cap.

    :param cap_type: The type of the cap.
    :type cap_type: CylinderCapType
    :param center_cap: The center of the cap.
    :type center_cap: Point3D
    :param radius_cap: The radius of the cap.
    :type radius_cap: float
    :param normal_cap: The normal vector of the cap.
    :type normal_cap: Vector3D
    :param center_cylinder: The center of the cylinder.
    :type center_cylinder: Point3D
    :param resolution: The resolution of the cap.
    :type resolution: int
    :param filter_for_pov: If True, only return points that are on the visible surface
        of the cap (i.e. the surface of the cap we can see from POV positive z-axis
        towards origin).
    :type filter_for_pov: bool
    :return: The points on the surface of the cap.
    :rtype: ty.List[Point3D]
    :raises ValueError: If the cap type is unknown.
    """
    if cap_type == CylinderCapType.NO_CAP:
        return []

    elif cap_type == CylinderCapType.FLAT:
        circle = Circle3D(center_cap, radius_cap, normal_cap)
        return get_points_on_circumference_circle_3d(circle, resolution)

    elif cap_type == CylinderCapType.ROUND:
        sphere = Sphere(center_cap, radius_cap)
        plane = Plane3D(center_cap, normal_cap)

        points = get_points_on_surface_sphere(
            sphere, resolution, resolution, filter_for_pov=filter_for_pov
        )

        points = [
            point for point in points if not same_side_of_plane(plane, center_cylinder, point)
        ]

        return points

    else:
        raise ValueError(f"Unknown cap type: '{cap_type}'")


def get_points_on_surface_cylinder(
    cylinder: Cylinder,
    resolution: int,
) -> ty.List[Point3D]:
    """Generate points on the surface of the cylinder.

    :param cylinder: The cylinder.
    :type cylinder: Cylinder
    :param resolution: The resolution of the cylinder.
    :type resolution: int
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

    cap_points = get_points_on_surface_cap(
        cap_type, cylinder.start, cylinder.radius, normal, cylinder.end, resolution, False
    )
    points.extend(cap_points)

    cap_points = get_points_on_surface_cap(
        cap_type, cylinder.end, cylinder.radius, normal, cylinder.start, resolution, False
    )
    points.extend(cap_points)

    return points


# ==============================================================================
# Check if points are inside a shape
# ==============================================================================


def point_is_inside_sphere(sphere: Sphere, point: Point3D) -> bool:
    """Check if a point is inside the sphere.

    :param sphere: The sphere to check.
    :type sphere: Sphere
    :param point: The point to check.
    :type point: Point3D
    :return: True if the point is inside the sphere, False otherwise.
    :rtype: bool
    """
    return sphere.center.calculate_distance(point) <= sphere.radius


def point_is_inside_cylinder(cylinder: Cylinder, point: Point3D) -> bool:
    """Check if a point is inside the cylinder.

    :param cylinder: The cylinder to check.
    :type cylinder: Cylinder
    :param point: The point to check.
    :type point: Point3D
    :return: True if the point is inside the cylinder, False otherwise.
    :rtype: bool
    :raises ValueError: If the cap type is unknown.
    """
    line = Line3D(cylinder.start, cylinder.end)
    dist = distance_to_line(line, point)
    cap_type = cylinder.cap_type

    if cap_type == CylinderCapType.ROUND:
        return dist <= cylinder.radius

    elif cap_type == CylinderCapType.FLAT or cap_type == CylinderCapType.NO_CAP:
        normal = cylinder.end.create_vector(cylinder.start).normalize()
        plane1 = Plane3D(cylinder.start, normal)
        plane2 = Plane3D(cylinder.end, normal)
        is_between_planes = same_side_of_plane(plane1, point, cylinder.end) and same_side_of_plane(
            plane2, point, cylinder.start
        )
        return dist <= cylinder.radius and is_between_planes

    else:
        raise ValueError(f"Unknown cap type: '{cap_type}'")


# ==============================================================================
# Check if shapes intersect
# ==============================================================================


def sphere_intersects_with_sphere(sphere1: Sphere, sphere2: Sphere) -> bool:
    """Check if two spheres intersect.

    :param sphere1: The first sphere.
    :type sphere1: Sphere
    :param sphere2: The second sphere.
    :type sphere2: Sphere
    :return: True if the spheres intersect, False otherwise.
    :rtype: bool
    """
    c1, r1 = sphere1.center, sphere1.radius
    c2, r2 = sphere2.center, sphere2.radius
    return c1.calculate_distance(c2) <= r1 + r2


def sphere_intersects_with_cylinder(sphere: Sphere, cylinder: Cylinder) -> bool:
    """Check if a sphere intersects with a cylinder.

    :param sphere: The sphere.
    :type sphere: Sphere
    :param cylinder: The cylinder.
    :type cylinder: Cylinder
    :return: True if the sphere intersects with the cylinder, False otherwise.
    :rtype: bool
    """
    d = sphere.radius + cylinder.radius
    return (
        sphere.center.calculate_distance(cylinder.start) <= d
        or sphere.center.calculate_distance(cylinder.end) <= d
    )


def cylinder_intersects_with_cylinder(cylinder1: Cylinder, cylinder2: Cylinder) -> bool:
    """Check if two cylinders intersect.

    :param cylinder1: The first cylinder.
    :type cylinder1: Cylinder
    :param cylinder2: The second cylinder.
    :type cylinder2: Cylinder
    :return: True if the cylinders intersect, False otherwise.
    :rtype: bool
    """
    d = cylinder1.radius + cylinder2.radius
    return (
        cylinder1.start.calculate_distance(cylinder2.start) <= d
        or cylinder1.start.calculate_distance(cylinder2.end) <= d
        or cylinder1.end.calculate_distance(cylinder2.start) <= d
        or cylinder1.end.calculate_distance(cylinder2.end) <= d
    )
