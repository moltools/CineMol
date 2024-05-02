# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.fitting module."""

import unittest

from cinemol.fitting import arange, argmax, argmin, calculate_convex_hull, process


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

