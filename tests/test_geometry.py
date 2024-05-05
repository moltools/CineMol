# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.geometry module."""

import math
import unittest

from cinemol.geometry import (
    Circle3D,
    Cylinder,
    CylinderCapType,
    Line3D,
    Plane3D,
    Point2D,
    Point3D,
    Sphere,
    Vector3D,
    cylinder_intersects_with_cylinder,
    distance_to_line,
    get_perpendicular_lines,
    get_points_on_circumference_circle_3d,
    get_points_on_line_3d,
    get_points_on_surface_cap,
    get_points_on_surface_circle_3d,
    get_points_on_surface_cylinder,
    get_points_on_surface_sphere,
    gram_schmidt,
    point_is_inside_cylinder,
    point_is_inside_sphere,
    same_side_of_plane,
    sign,
    sphere_intersects_with_cylinder,
    sphere_intersects_with_sphere,
)


class TestVector3D(unittest.TestCase):
    """Test the Vector3D class."""

    def test_vector3d_creation(self):
        """Test the creation of a Vector3D object."""
        vector = Vector3D(1, 2, 3)
        self.assertEqual(vector.x, 1)
        self.assertEqual(vector.y, 2)
        self.assertEqual(vector.z, 3)

    def test_random_vector3d_creation(self):
        """Test the creation of a Vector3D object with random values."""
        vector = Vector3D.create_random()
        self.assertIsNotNone(vector)

        # All values should be between 0 and -1.
        self.assertTrue(0 <= vector.x <= 1)
        self.assertTrue(0 <= vector.y <= 1)
        self.assertTrue(0 <= vector.z <= 1)

    def test_vector3d_length(self):
        """Test the length of a Vector3D object."""
        vector = Vector3D(1, 2, 3)
        self.assertAlmostEqual(vector.length(), 3.742, places=3)

    def test_vector3d_normalize(self):
        """Test the normalization of a Vector3D object."""
        vector = Vector3D(1, 2, 3)
        normalized_vector = vector.normalize()
        self.assertAlmostEqual(normalized_vector.length(), 1.0, places=3)

    def test_vector3d_dot_product(self):
        """Test the dot product of two Vector3D objects."""
        vector1 = Vector3D(1, 2, 3)
        vector2 = Vector3D(4, 5, 6)
        result = vector1.dot(vector2)
        self.assertEqual(result, 32)

    def test_vector3d_cross_product(self):
        """Test the cross product of two Vector3D objects."""
        vector1 = Vector3D(1, 2, 3)
        vector2 = Vector3D(4, 5, 6)
        result = vector1.cross(vector2)
        self.assertEqual(result.x, -3)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.z, -3)

    def test_vector3d_subtract(self):
        """Test the subtraction of two Vector3D objects."""
        vector1 = Vector3D(1, 2, 3)
        vector2 = Vector3D(4, 5, 6)
        result = vector1.subtract(vector2)
        self.assertEqual(result.x, -3)
        self.assertEqual(result.y, -3)
        self.assertEqual(result.z, -3)

    def test_vector3d_multiply(self):
        """Test the multiplication of a Vector3D object by a scalar."""
        vector = Vector3D(1, 2, 3)
        result = vector.multiply(2)
        self.assertEqual(result.x, 2)
        self.assertEqual(result.y, 4)
        self.assertEqual(result.z, 6)


class TestPoint2D(unittest.TestCase):
    """Test the Point2D class."""

    def test_point2d_creation(self):
        """Test the creation of a Point2D object."""
        point = Point2D(1, 2)
        self.assertEqual(point.x, 1)
        self.assertEqual(point.y, 2)

    def test_poin2d_subtract(self):
        """Test the subtraction of two Point2D objects."""
        point1 = Point2D(1, 2)
        point2 = Point2D(3, 4)
        result = point1.subtract_point(point2)
        self.assertEqual(result.x, -2)
        self.assertEqual(result.y, -2)

    def test_point2d_cross_product(self):
        """Test the cross product of two Point2D objects."""
        point1 = Point2D(1, 2)
        point2 = Point2D(3, 4)
        result = point1.cross(point2)
        self.assertEqual(result, -2)


class TestPoint3D(unittest.TestCase):
    """Test the Point3D class."""

    def test_point3d_creation(self):
        """Test the creation of a Point3D object."""
        point = Point3D(1, 2, 3)
        self.assertEqual(point.x, 1)
        self.assertEqual(point.y, 2)
        self.assertEqual(point.z, 3)

    def test_create_vector(self):
        """Test the creation of a Vector3D object from two Point3D objects."""
        point1 = Point3D(1, 2, 3)
        point2 = Point3D(4, 5, 6)
        vector = point1.create_vector(point2)
        self.assertEqual(vector.x, 3)
        self.assertEqual(vector.y, 3)
        self.assertEqual(vector.z, 3)

    def test_calculate_distance(self):
        """Test the calculation of the distance between two Point3D objects."""
        point1 = Point3D(1, 2, 3)
        point2 = Point3D(4, 5, 6)
        distance = point1.calculate_distance(point2)
        self.assertAlmostEqual(distance, 5.196, places=3)

    def test_calculate_midpoint(self):
        """Test the calculation of the midpoint between two Point3D objects."""
        point1 = Point3D(1, 2, 3)
        point2 = Point3D(4, 5, 6)
        midpoint = point1.midpoint(point2)
        self.assertEqual(midpoint.x, 2.5)
        self.assertEqual(midpoint.y, 3.5)
        self.assertEqual(midpoint.z, 4.5)

    def test_translate_with_vector(self):
        """Test the translation of a Point3D object with a Vector3D object."""
        point = Point3D(1, 2, 3)
        vector = Vector3D(4, 5, 6)
        result = point.translate(vector)
        self.assertEqual(result.x, 5)
        self.assertEqual(result.y, 7)
        self.assertEqual(result.z, 9)

    def test_rotate_around_x_axis(self):
        """Test the rotation of a Point3D object around an axis, in radians."""
        point = Point3D(1, 0, 0)
        rotation = 90 * math.pi / 180
        rotated_point = point.rotate(rotation)

        # The point is on the x-axis, so the position remains the same.
        self.assertAlmostEqual(rotated_point.x, 1, places=3)
        self.assertAlmostEqual(rotated_point.y, 0, places=3)
        self.assertAlmostEqual(rotated_point.z, 0, places=3)

    def test_rotate_around_x_axis_with_offset(self):
        """Test the rotation of a Point3D object around an axis, in radians, with an offset."""
        point = Point3D(1, 1, 0)
        rotation = 90 * math.pi / 180
        rotated_point = point.rotate(rotation)

        # The point is on the x-axis, so the position remains the same.
        self.assertAlmostEqual(rotated_point.x, 1, places=3)
        self.assertAlmostEqual(rotated_point.y, 0, places=3)
        self.assertAlmostEqual(rotated_point.z, 1, places=3)


class TestSign(unittest.TestCase):
    """Test the sign function."""

    def test_sign_positive(self):
        """Test the sign function with a positive value."""
        result = sign(1)
        self.assertEqual(result, 1)

    def test_sign_negative(self):
        """Test the sign function with a negative value."""
        result = sign(-1)
        self.assertEqual(result, -1)

    def test_sign_zero(self):
        """Test the sign function with a zero value."""
        result = sign(0)
        self.assertEqual(result, 0)


class TestGramSchmidt(unittest.TestCase):
    """Test the Gram-Schmidt process."""

    def test_gram_schmidt(self):
        """Test the Gram-Schmidt process."""
        vector_n = Vector3D(1, 2, 3)
        result = gram_schmidt(vector_n)
        self.assertIsNotNone(result)
        vector_v, vector_w = result
        self.assertAlmostEqual(vector_n.dot(vector_w), 0, places=3)
        self.assertAlmostEqual(vector_v.dot(vector_w), 0, places=3)


class TestLine3D(unittest.TestCase):
    """Test the Line3D class."""

    def test_line3d_creation(self):
        """Test the creation of a Line3D object."""
        point1 = Point3D(1, 2, 3)
        point2 = Point3D(4, 5, 6)
        line = Line3D(point1, point2)
        self.assertEqual(line.start, point1)
        self.assertEqual(line.end, point2)


class TestPlane3D(unittest.TestCase):
    """Test the Plane3D class."""

    def test_plane3d_creation(self):
        """Test the creation of a Plane3D object."""
        center = Point3D(1, 2, 3)
        normal = Vector3D(4, 5, 6)
        plane = Plane3D(center, normal)
        self.assertEqual(plane.center, center)
        self.assertEqual(plane.normal, normal)


class TestCircle3D(unittest.TestCase):
    """Test the Circle3D class."""

    def test_circle3d_creation(self):
        """Test the creation of a Circle3D object."""
        center = Point3D(1, 2, 3)
        normal = Vector3D(4, 5, 6)
        radius = 7
        circle = Circle3D(center, radius, normal)
        self.assertEqual(circle.center, center)
        self.assertEqual(circle.normal, normal)
        self.assertEqual(circle.radius, radius)


class TestSphere(unittest.TestCase):
    """Test the Sphere class."""

    def test_sphere_creation(self):
        """Test the creation of a Sphere object."""
        center = Point3D(1, 2, 3)
        radius = 4
        sphere = Sphere(center, radius)
        self.assertEqual(sphere.center, center)
        self.assertEqual(sphere.radius, radius)


class TestCylinderCapType(unittest.TestCase):
    """Test the CylinderCapType class."""

    def test_cylinder_cap_type_creation(self):
        """Test the creation of a CylinderCapType object."""
        cap_type = CylinderCapType.NO_CAP
        self.assertEqual(cap_type, CylinderCapType.NO_CAP)
        cap_type = CylinderCapType.FLAT
        self.assertEqual(cap_type, CylinderCapType.FLAT)
        cap_type = CylinderCapType.ROUND
        self.assertEqual(cap_type, CylinderCapType.ROUND)


class TestCylinder(unittest.TestCase):
    """Test the Cylinder class."""

    def test_cylinder_creation(self):
        """Test the creation of a Cylinder object."""
        start = Point3D(1, 2, 3)
        end = Point3D(4, 5, 6)
        radius = 7
        cap_type = CylinderCapType.NO_CAP
        cylinder = Cylinder(start, end, radius, cap_type)
        self.assertEqual(cylinder.start, start)
        self.assertEqual(cylinder.end, end)
        self.assertEqual(cylinder.radius, radius)
        self.assertEqual(cylinder.cap_type, cap_type)


class TestSameSideOfPlane(unittest.TestCase):
    """Test the same_side_of_plane function."""

    def test_same_side_of_plane_same_side(self):
        """Test same_side_of_plane with points on the same side."""
        plane_point = Point3D(0, 0, 0)
        plane_normal = Vector3D(0, 0, 1)
        plane = Plane3D(plane_point, plane_normal)
        point1 = Point3D(1, 1, 1)
        point2 = Point3D(-1, -1, 1)
        result = same_side_of_plane(plane, point1, point2)
        self.assertTrue(result)

    def test_same_side_of_plane_opposite_side(self):
        """Test same_side_of_plane with points on opposite sides."""
        plane_point = Point3D(0, 0, 0)
        plane_normal = Vector3D(0, 0, 1)
        plane = Plane3D(plane_point, plane_normal)
        point1 = Point3D(1, 1, 1)
        point2 = Point3D(-1, -1, -1)
        result = same_side_of_plane(plane, point1, point2)
        self.assertFalse(result)


class TestDistanceToLine(unittest.TestCase):
    """Test the distance_to_line function."""

    def test_distance_to_line_is_zero(self):
        """Test the distance_to_line function."""
        point = Point3D(1, 1, 1)
        line = Line3D(Point3D(0, 0, 0), Point3D(2, 2, 2))
        result = distance_to_line(line, point)
        self.assertAlmostEqual(result, 0.0, places=3)

    def test_distance_to_line_is_not_zero(self):
        """Test the distance_to_line function."""
        point = Point3D(1, 0, 0)
        line = Line3D(Point3D(-1, -1, 0), Point3D(2, 2, 0))
        result = distance_to_line(line, point)
        self.assertAlmostEqual(result, 0.707, places=3)


class TestGetPerpendicularLines(unittest.TestCase):
    """Test the get_perpendicular_lines function."""

    def test_zero_perpendicular_lines(self):
        """Returning zero lines raises an exception."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        width = 1
        with self.assertRaises(ValueError):
            get_perpendicular_lines(line, width, 0)

    def test_negative_perpendicular_lines(self):
        """Returning negative lines raises an exception."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        width = 1
        with self.assertRaises(ValueError):
            get_perpendicular_lines(line, width, -1)

    def test_zero_width(self):
        """Returning zero width raises an exception."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        width = 0
        with self.assertRaises(ValueError):
            get_perpendicular_lines(line, width, 1)

    def test_negative_width(self):
        """Returning negative width raises an exception."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        width = -1
        with self.assertRaises(ValueError):
            get_perpendicular_lines(line, width, 1)

    def test_one_perpendicular_line(self):
        """Test one perpendicular line."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        width = 1
        result = get_perpendicular_lines(line, width, 1)
        self.assertEqual([line], result)

    def test_two_perpendicular_lines(self):
        """Test two perpendicular lines."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0))
        width = 1
        result = get_perpendicular_lines(line, width, 2)

        # Distances between the starts and the ends of the two new lines should
        # both be equal to the given width.
        self.assertAlmostEqual(result[0].start.calculate_distance(result[1].start), width, places=3)
        self.assertAlmostEqual(result[0].end.calculate_distance(result[1].end), width, places=3)

        # Distances between the starts and the ends of the two new lines and the
        # original line should be equal to half the given width.
        self.assertAlmostEqual(result[0].start.calculate_distance(line.start), width / 2, places=3)
        self.assertAlmostEqual(result[0].end.calculate_distance(line.end), width / 2, places=3)

    def test_three_perpendicular_lines(self):
        """Test three perpendicular lines."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0))
        width = 1
        result = get_perpendicular_lines(line, width, 3)

        # Distances between the starts and the ends of the three new lines should
        # all be equal to the given width.
        self.assertAlmostEqual(result[0].start.calculate_distance(result[1].start), width, places=3)
        self.assertAlmostEqual(result[0].end.calculate_distance(result[1].end), width, places=3)
        self.assertAlmostEqual(result[1].start.calculate_distance(result[2].start), width, places=3)
        self.assertAlmostEqual(result[1].end.calculate_distance(result[2].end), width, places=3)
        self.assertAlmostEqual(
            result[0].start.calculate_distance(result[2].start), width * 2, places=3
        )
        self.assertAlmostEqual(result[0].end.calculate_distance(result[2].end), width * 2, places=3)

        # Distances between the starts and the ends of the three new lines and the
        # original line should be equal to the given width for the first and third line ...
        self.assertAlmostEqual(result[0].start.calculate_distance(line.start), width, places=3)
        self.assertAlmostEqual(result[0].end.calculate_distance(line.end), width, places=3)
        self.assertAlmostEqual(result[2].start.calculate_distance(line.start), width, places=3)
        self.assertAlmostEqual(result[2].end.calculate_distance(line.end), width, places=3)

        # ... and zero for the second line.
        self.assertAlmostEqual(result[1].start.calculate_distance(line.start), 0, places=3)
        self.assertAlmostEqual(result[1].end.calculate_distance(line.end), 0, places=3)


class TestGetPointsOnLine3D(unittest.TestCase):
    """Test the get_points_on_line_3d function."""

    def test_get_exact_numer_of_points_with_get_points_on_line_3d(self):
        """Test if the get_points_on_line_3d function returns the correct number of points."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        result = get_points_on_line_3d(line, 10)
        self.assertEqual(len(result), 11)
        self.assertAlmostEqual(result[0].calculate_distance(line.start), 0, places=3)
        self.assertAlmostEqual(result[-1].calculate_distance(line.end), 0, places=3)

    def test_points_on_line_3d_returns_points_on_line(self):
        """Test if the get_points_on_line_3d function returns points on the line."""
        line = Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        result = get_points_on_line_3d(line, 10)
        for i in range(1, len(result)):
            self.assertAlmostEqual(distance_to_line(line, result[i]), 0, places=3)


class TestGetPointsOnCircumferenceCircle3D(unittest.TestCase):
    """Test the get_points_on_circumference_circle_3d function."""

    def test_get_exact_number_of_points_with_get_points_on_circumference_circle_3d(self):
        """Test if the get_points_on_circumference_circle_3d function returns the correct number of points."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        result = get_points_on_circumference_circle_3d(circle, 10)
        self.assertEqual(len(result), 10)

    def test_points_on_circumference_circle_3d_returns_points_on_circle(self):
        """Test if the get_points_on_circumference_circle_3d function returns points on the circle."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        result = get_points_on_circumference_circle_3d(circle, 10)
        for point in result:
            self.assertAlmostEqual(point.calculate_distance(circle.center), circle.radius, places=3)


class TestGetPointsOnSurfaceCircle3D(unittest.TestCase):
    """Test the get_points_on_surface_circle_3d function."""

    def test_get_exact_number_of_points_with_get_points_on_surface_circle_3d(self):
        """Test if the get_points_on_surface_circle_3d function returns the correct number of points."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        num_radii = 3
        num_points = 10
        result = get_points_on_surface_circle_3d(circle, num_radii, num_points)
        self.assertEqual(len(result), num_radii * num_points)

    def test_points_on_surface_circle_3d_returns_points_on_surface(self):
        """Test if the get_points_on_surface_circle_3d function returns points on the surface."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        num_radii = 3
        num_points = 10
        result = get_points_on_surface_circle_3d(circle, num_radii, num_points)
        for point in result:
            self.assertTrue(point.calculate_distance(circle.center) <= circle.radius)
            self.assertAlmostEqual(point.z, 0, places=3)


class TestGetPointsOnSurfaceSphere(unittest.TestCase):
    """Test the get_points_on_surface_sphere function."""

    def test_points_on_surface_sphere_returns_points_on_surface(self):
        """Test if the get_points_on_surface_sphere function returns points on the surface."""
        sphere = Sphere(Point3D(0, 0, 0), 1)
        result = get_points_on_surface_sphere(sphere, 10, 10)
        for point in result:
            self.assertAlmostEqual(point.calculate_distance(sphere.center), sphere.radius, places=3)


class TestGetPointsOnSurfaceCap(unittest.TestCase):
    """Test the get_points_on_surface_cap function."""

    def test_points_on_surface_flat_cap_returns_points_on_surface(self):
        """Test if the get_points_on_surface_cap function returns points on the surface."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        cap_type = CylinderCapType.FLAT
        result = get_points_on_surface_cap(
            cap_type, circle.center, circle.radius, circle.normal, Point3D(0, 0, -1), 10
        )
        for point in result:
            self.assertTrue(round(point.calculate_distance(circle.center), 5) <= circle.radius)
            self.assertAlmostEqual(point.z, 0, places=3)

    def test_points_on_surface_round_cap_returns_points_on_surface(self):
        """Test if the get_points_on_surface_cap function returns points on the surface."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        cap_type = CylinderCapType.ROUND
        result = get_points_on_surface_cap(
            cap_type, circle.center, circle.radius, circle.normal, Point3D(0, 0, -1), 10
        )
        for point in result:
            self.assertAlmostEqual(point.calculate_distance(circle.center), circle.radius, places=3)
            self.assertTrue(point.z >= 0.0)

    def test_points_on_surface_no_cap_returns_no_points(self):
        """Test if the get_points_on_surface_cap function returns no points."""
        circle = Circle3D(Point3D(0, 0, 0), 1, Vector3D(0, 0, 1))
        cap_type = CylinderCapType.NO_CAP
        result = get_points_on_surface_cap(
            cap_type, circle.center, circle.radius, circle.normal, Point3D(0, 0, -1), 10
        )
        self.assertEqual(len(result), 0)


class TestGetPointsOnSurfaceCylinder(unittest.TestCase):
    """Test the get_points_on_surface_cylinder function."""

    def test_points_on_surface_cylinder_returns_points_on_surface(self):
        """Test if the get_points_on_surface_cylinder function returns points on the surface."""
        start = Point3D(0, 0, 0)
        end = Point3D(0, 0, 1)
        line = Line3D(start, end)
        radius = 1
        cylinder = Cylinder(start, end, radius, CylinderCapType.NO_CAP)
        result = get_points_on_surface_cylinder(cylinder, 10)
        for point in result:
            self.assertAlmostEqual(distance_to_line(line, point), radius, places=3)


class TestPointIsInsideSphere(unittest.TestCase):
    """Test the point_is_inside_sphere function."""

    def test_point_is_inside_sphere_inside(self):
        """Test if the point_is_inside_sphere function returns True for a point inside the sphere."""
        sphere = Sphere(Point3D(0, 0, 0), 1)
        point = Point3D(0, 0, 0.5)
        result = point_is_inside_sphere(sphere, point)
        self.assertTrue(result)

    def test_point_is_inside_sphere_outside(self):
        """Test if the point_is_inside_sphere function returns False for a point outside the sphere."""
        sphere = Sphere(Point3D(0, 0, 0), 1)
        point = Point3D(0, 0, 1.5)
        result = point_is_inside_sphere(sphere, point)
        self.assertFalse(result)


class TestPointIsInsideCylinder(unittest.TestCase):
    """Test the point_is_inside_cylinder function."""

    def test_point_is_inside_cylinder_inside(self):
        """Test if the point_is_inside_cylinder function returns True for a point inside the cylinder."""
        cylinder = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1, CylinderCapType.NO_CAP)
        point = Point3D(0, 0, 0.5)
        result = point_is_inside_cylinder(cylinder, point)
        self.assertTrue(result)

    def test_point_is_inside_cylinder_outside(self):
        """Test if the point_is_inside_cylinder function returns False for a point outside the cylinder."""
        cylinder = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1, CylinderCapType.NO_CAP)
        point = Point3D(0, 0, 1.5)
        result = point_is_inside_cylinder(cylinder, point)
        self.assertFalse(result)


class TestSphereIntersectsWithSphere(unittest.TestCase):
    """Test the sphere_intersects_with_sphere function."""

    def test_sphere_intersects_with_sphere_intersects(self):
        """Test if the sphere_intersects_with_sphere function returns True for intersecting spheres."""
        sphere1 = Sphere(Point3D(0, 0, 0), 1)
        sphere2 = Sphere(Point3D(1, 0, 0), 1)
        result = sphere_intersects_with_sphere(sphere1, sphere2)
        self.assertTrue(result)

    def test_sphere_intersects_with_sphere_does_not_intersect(self):
        """Test if the sphere_intersects_with_sphere function returns False for non-intersecting spheres."""
        sphere1 = Sphere(Point3D(0, 0, 0), 1)
        sphere2 = Sphere(Point3D(3, 0, 0), 1)
        result = sphere_intersects_with_sphere(sphere1, sphere2)
        self.assertFalse(result)


class TestSphereIntersectsWithCylinder(unittest.TestCase):
    """Test the sphere_intersects_with_cylinder function."""

    def test_sphere_intersects_with_cylinder_intersects(self):
        """Test if sphere_intersects_with_cylinder function returns True for intersecting sphere and cylinder."""
        sphere = Sphere(Point3D(0, 0, 0), 1)
        cylinder = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1, CylinderCapType.NO_CAP)
        result = sphere_intersects_with_cylinder(sphere, cylinder)
        self.assertTrue(result)

    def test_sphere_intersects_with_cylinder_does_not_intersect(self):
        """Test if sphere_intersects_with_cylinder function returns False for non-intersecting sphere and cylinder."""
        sphere = Sphere(Point3D(0, 0, 0), 1)
        cylinder = Cylinder(Point3D(0, 0, 3), Point3D(0, 0, 4), 1, CylinderCapType.NO_CAP)
        result = sphere_intersects_with_cylinder(sphere, cylinder)
        self.assertFalse(result)


class TestCylinderIntersectsWithCylinder(unittest.TestCase):
    """Test the cylinder_intersects_with_cylinder function."""

    def test_cylinder_intersects_with_cylinder_intersects(self):
        """Test if cylinder_intersects_with_cylinder function returns True for intersecting cylinders."""
        cylinder1 = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1, CylinderCapType.NO_CAP)
        cylinder2 = Cylinder(Point3D(1, 0, 0), Point3D(1, 0, 1), 1, CylinderCapType.NO_CAP)
        result = cylinder_intersects_with_cylinder(cylinder1, cylinder2)
        self.assertTrue(result)

    def test_cylinder_intersects_with_cylinder_does_not_intersect(self):
        """Test if cylinder_intersects_with_cylinder function returns False for non-intersecting cylinders."""
        cylinder1 = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1, CylinderCapType.NO_CAP)
        cylinder2 = Cylinder(Point3D(3, 0, 0), Point3D(3, 0, 1), 1, CylinderCapType.NO_CAP)
        result = cylinder_intersects_with_cylinder(cylinder1, cylinder2)
        self.assertFalse(result)
