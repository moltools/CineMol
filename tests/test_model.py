# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.model module."""

import unittest

from cinemol.geometry import CylinderCapType, Cylinder, Line3D, Sphere, Point3D
from cinemol.model import (
    ModelCylinder,
    ModelSphere,
    ModelWire,
    Scene,
    calculate_intersecting_nodes,
    create_fill,
    get_node_polygon_vertices,
    prepare_nodes_for_intersecting,
)
from cinemol.style import Cartoon, Color


class TestModelSphere(unittest.TestCase):
    """Test the ModelSphere class."""

    def test_init(self):
        """Test ModelSphere creation."""
        sphere = Sphere(Point3D(1, 2, 3), 1.0)
        depiction = Cartoon(Color(0, 255, 0))
        model = ModelSphere(sphere, depiction)
        self.assertEqual(model.geometry, sphere)
        self.assertEqual(model.depiction, depiction)


class TestModelCylinder(unittest.TestCase):
    """Test the ModelCylinder class."""

    def test_init(self):
        """Test ModelCylinder creation."""
        cylinder = Cylinder(
            Point3D(0, 0, 0),
            Point3D(0, 0, 1),
            1.0,
            CylinderCapType.NO_CAP
        )
        depiction = Cartoon(Color(0, 255, 0))
        model = ModelCylinder(cylinder, depiction)
        self.assertEqual(model.geometry, cylinder)
        self.assertEqual(model.depiction, depiction)


class TestModelWire(unittest.TestCase):
    """Test the ModelWire class."""

    def test_init(self):
        """Test ModelWire creation."""
        line = Line3D(Point3D(0, 0, 0), Point3D(0, 0, 1))
        color = Color(0, 255, 0)
        width = 1.0
        opacity = 0.5
        model = ModelWire(line, color, width, opacity)
        self.assertEqual(model.geometry, line)
        self.assertEqual(model.color, color)
        self.assertEqual(model.width, width)
        self.assertEqual(model.opacity, opacity)
