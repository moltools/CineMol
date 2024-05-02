# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.fitting module."""

import unittest

from cinemol.fitting import arange, argmax, argmin, calculate_convex_hull, process
from cinemol.geometry import Point2D


class TestArgmin(unittest.TestCase):
    """Test the argmin function."""

    def test_argmin_basic(self):
        """Test the argmin with a basic input."""
        result = argmin([4, 2, 9, 3])
        self.assertEqual(result, 1)

    def test_argmin_all_same(self):
        """Test argmin when all elements are the same."""
        result = argmin([5, 5, 5, 5])
        self.assertEqual(result, 0)

    def test_argmin_negative(self):
        """Test argmin with only negative numbers."""
        result = argmin([-2, -3, -1, -4])
        self.assertEqual(result, 3)

    def test_argmin_mixed_positive_and_negative(self):
        """Test argmin with mixed positive and negative numbers."""
        result = argmin([3, -1, 4, -2])
        self.assertEqual(result, 3)


class TestArgmax(unittest.TestCase):
    """Test the argmax function."""

    def test_argmax_basic(self):
        """Test the argmax with a basic input."""
        result = argmax([4, 2, 9, 3])
        self.assertEqual(result, 2)

    def test_argmax_all_same(self):
        """Test argmax when all elements are the same."""
        result = argmax([5, 5, 5, 5])
        self.assertEqual(result, 0)

    def test_argmax_negative(self):
        """Test argmax with only negative numbers."""
        result = argmax([-2, -3, -1, -4])
        self.assertEqual(result, 2)

    def test_argmax_mixed_positive_and_negative(self):
        """Test argmax with mixed positive and negative numbers."""
        result = argmax([3, -1, 4, -2])
        self.assertEqual(result, 2)


class TestArange(unittest.TestCase):
    """Test the arange function."""

    def test_arange_basic(self):
        """Test the arange with a basic input."""
        result = arange(5)
        self.assertEqual(result, [0, 1, 2, 3, 4])

    def test_arange_negative(self):
        """Test the arange with a negative input."""
        result = arange(-3)
        self.assertEqual(result, [])

    def test_arange_zero(self):
        """Test the arange with a zero input."""
        result = arange(0)
        self.assertEqual(result, [])

    def test_arange_float(self):
        """Test the arange with a float input. This should raise a TypeError."""
        with self.assertRaises(TypeError):
            arange(3.5)

    def test_arange_negative_float(self):
        """Test the arange with a negative float input. This should raise a TypeError."""
        with self.assertRaises(TypeError):
            arange(-3.5)


class TestProcess(unittest.TestCase):
    """Test the process function."""

    def test_process_all_points_on_positive_side(self):
        """Test the process for points on the positive side of the line."""
        points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        indices_of_points = arange(len(points))
        index_a = 0
        index_b = 2
        result = process(points, indices_of_points, index_a, index_b)
        self.assertEqual(result, [0, 3, 2])

    def test_process_all_points_on_negative_side(self):
        """Test the process for points on the negative side of the line."""
        points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        indices_of_points = arange(len(points))
        index_a = 0
        index_b = 2
        result = process(points, indices_of_points, index_b, index_a)
        self.assertEqual(result, [2, 1, 0])

    def test_process_mixed_points(self):
        """Test the process with mixed points on the positive and negative sides of the line."""
        points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        indices_of_points = arange(len(points))
        index_a = 0
        index_b = 2
        result = process(points, indices_of_points, index_a, index_b)
        self.assertEqual(result, [0, 3, 2])

    def test_process_with_incorrect_points_input(self):
        """Test the process with incorrect points input. This should raise an AttributeError."""
        points = [(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        indices_of_points = arange(len(points))
        index_a = 0
        index_b = 2
        with self.assertRaises(AttributeError):
            process(points, indices_of_points, index_b, index_a)


class TestCalculateConvexHull(unittest.TestCase):
    """Test the calculate_convex_hull function."""

    def test_calculate_convex_hull_basic(self):
        """Test the calculate_convex_hull with a basic input."""
        points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        result = calculate_convex_hull(points)
        self.assertEqual(result, [0, 3, 2, 1])

    def test_calculate_convex_hull_including_centroid(self):
        """Test the calculate_convex_hull with the centroid included."""
        points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1), Point2D(1, 0)]
        result = calculate_convex_hull(points)
        self.assertEqual(result, [0, 3, 2, 1])

    def test_calculate_convex_hull_incorrect_points_input(self):
        """Test the calculate_convex_hull with incorrect points input. This should raise an AttributeError."""
        points = [(0, 0), Point2D(1, 1), Point2D(2, 0), Point2D(1, -1)]
        with self.assertRaises(AttributeError):
            calculate_convex_hull(points)
