# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.model module."""

import unittest

from cinemol.geometry import Cylinder, CylinderCapType, Line3D, Point3D, Sphere
from cinemol.model import ModelCylinder, ModelSphere, ModelWire, Scene
from cinemol.style import Cartoon, Color


class TestModelSphere(unittest.TestCase):
    """Test the ModelSphere class."""

    def test_creation_model_sphere(self):
        """Test ModelSphere creation."""
        sphere = Sphere(Point3D(1, 2, 3), 1.0)
        depiction = Cartoon(Color(0, 255, 0))
        model = ModelSphere(sphere, depiction)
        self.assertEqual(model.geometry, sphere)
        self.assertEqual(model.depiction, depiction)


class TestModelCylinder(unittest.TestCase):
    """Test the ModelCylinder class."""

    def test_creation_model_cylinder(self):
        """Test ModelCylinder creation."""
        cylinder = Cylinder(Point3D(0, 0, 0), Point3D(0, 0, 1), 1.0, CylinderCapType.NO_CAP)
        depiction = Cartoon(Color(0, 255, 0))
        model = ModelCylinder(cylinder, depiction)
        self.assertEqual(model.geometry, cylinder)
        self.assertEqual(model.depiction, depiction)


class TestModelWire(unittest.TestCase):
    """Test the ModelWire class."""

    def test_creation_model_wire(self):
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


class TestScene(unittest.TestCase):
    """Test the Scene class."""

    def test_creation_empty_scene(self):
        """Test creation of Scene without supplying nodes."""
        scene = Scene()
        self.assertEqual(scene.nodes, [])

    def test_createion_scene_with_single_node(self):
        """Test creation of Scene with a single node."""
        sphere = Sphere(Point3D(1, 2, 3), 1.0)
        depiction = Cartoon(Color(0, 255, 0))
        model = ModelSphere(sphere, depiction)
        scene = Scene([model])
        self.assertEqual(scene.nodes, [model])
