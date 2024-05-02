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
    gram_schmidt,
    sign,
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
