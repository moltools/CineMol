# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.api module."""

import unittest

from cinemol.api import (
    Atom,
    Bond,
    Look,
    Style,
    draw_atoms_in_spacefilling_style,
    draw_bonds_in_tube_style,
    draw_bonds_in_wireframe_style,
    draw_molecule,
    filter_atoms_and_bonds,
)
from cinemol.geometry import CylinderCapType
from cinemol.model import Scene
from cinemol.style import Color


class TestStyle(unittest.TestCase):
    """Test the Style class."""

    def test_style(self):
        """Test the creation of a Style object."""
        style_type = Style.SPACEFILLING
        self.assertEqual(style_type, Style.SPACEFILLING)
        style_type = Style.BALL_AND_STICK
        self.assertEqual(style_type, Style.BALL_AND_STICK)
        style_type = Style.TUBE
        self.assertEqual(style_type, Style.TUBE)
        style_type = Style.WIREFRAME
        self.assertEqual(style_type, Style.WIREFRAME)


class TestLook(unittest.TestCase):
    """Test the Look class."""

    def test_look(self):
        """Test the creation of a Look object."""
        look_type = Look.CARTOON
        self.assertEqual(look_type, Look.CARTOON)
        look_type = Look.GLOSSY
        self.assertEqual(look_type, Look.GLOSSY)


class TestAtom(unittest.TestCase):
    """Test the Atom class."""

    def test_creation_carbon_atom_with_defaults(self):
        """Test the creation of a carbon type Atom object with defaults."""
        atom = Atom(0, "C", (0, 0, 0))
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.symbol, "C")
        self.assertEqual(atom.coordinates, (0, 0, 0))
        self.assertEqual(atom.radius, None)
        self.assertEqual(atom.color, None)
        self.assertEqual(atom.opacity, 1.0)

    def test_creation_carbon_atom_without_defaults(self):
        """Test the creation of a carbon type Atom object without defaults."""
        atom = Atom(0, "C", (0, 0, 0), 1.0, (255, 0, 0), 0.5)
        self.assertEqual(atom.index, 0)
        self.assertEqual(atom.symbol, "C")
        self.assertEqual(atom.coordinates, (0, 0, 0))
        self.assertEqual(atom.radius, 1.0)
        self.assertEqual(atom.color, (255, 0, 0))
        self.assertEqual(atom.opacity, 0.5)


class TestBond(unittest.TestCase):
    """Test the Bond class."""

    def test_creation_single_bond_with_defaults(self):
        """Test the creation of a single type Bond object with defaults."""
        bond = Bond(0, 1, 1)
        self.assertEqual(bond.start_index, 0)
        self.assertEqual(bond.end_index, 1)
        self.assertEqual(bond.order, 1)
        self.assertEqual(bond.radius, None)
        self.assertEqual(bond.color, None)
        self.assertEqual(bond.opacity, 1.0)

    def test_creation_single_bond_without_defaults(self):
        """Test the creation of a single type Bond object without defaults."""
        bond = Bond(0, 1, 1, 0.5, (255, 0, 0), 0.5)
        self.assertEqual(bond.start_index, 0)
        self.assertEqual(bond.end_index, 1)
        self.assertEqual(bond.order, 1)
        self.assertEqual(bond.radius, 0.5)
        self.assertEqual(bond.color, (255, 0, 0))
        self.assertEqual(bond.opacity, 0.5)


class TestFilterAtomsAndBonds(unittest.TestCase):
    """Test the filter_atoms_and_bonds function."""

    def test_filter_atoms_and_bonds(self):
        """Test the filter_atoms_and_bonds function without filtering."""
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        filtered_atoms, filtered_bonds = filter_atoms_and_bonds(atoms, bonds)
        self.assertEqual(len(filtered_atoms), 4)
        self.assertEqual(len(filtered_bonds), 3)

    def test_filter_atoms_and_bonds_with_filtering(self):
        """Test the filter_atoms_and_bonds function with filtering."""
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        filtered_atoms, filtered_bonds = filter_atoms_and_bonds(atoms, bonds, ["H"])
        self.assertEqual(len(filtered_atoms), 1)
        self.assertEqual(len(filtered_bonds), 0)

    def test_filter_atoms_and_bonds_with_filtering_out_all_atoms(self):
        """Test the filter_atoms_and_bonds function with filtering out multiple atom types."""
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        filtered_atoms, filtered_bonds = filter_atoms_and_bonds(atoms, bonds, ["C", "H"])
        self.assertEqual(len(filtered_atoms), 0)
        self.assertEqual(len(filtered_bonds), 0)


class TestDrawBondsInWireframeStyle(unittest.TestCase):
    """Test the draw_bonds_in_wireframe_style function."""

    def test_draw_bonds_in_wireframe_style_without_supplying_bonds(self):
        """Draw simple molecule in wireframe style without supplying bonds."""
        scene = Scene()
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        draw_bonds_in_wireframe_style(scene, atoms, [])
        self.assertEqual(len(scene.nodes), 0)

    def test_draw_bonds_in_wireframe_style_without_supplying_atoms(self):
        """Draw simple molecule in wireframe style without supplying atoms throws KeyError."""
        scene = Scene()
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        with self.assertRaises(KeyError):
            draw_bonds_in_wireframe_style(scene, [], bonds)

    def test_draw_bonds_in_wireframe_style(self):
        """Draw simple molecule in wireframe style."""
        scene = Scene()
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        draw_bonds_in_wireframe_style(scene, atoms, bonds)
        self.assertEqual(len(scene.nodes), 6)  # 3 x 2 half bonds.


class TtestDrawAtomsInSpacefillingStyle(unittest.TestCase):
    """Test the draw_atoms_in_spacefilling_style function."""

    def test_draw_atoms_in_spacefilling_style_without_supplying_atoms(self):
        """Draw simple molecule in spacefilling style without supplying atoms."""
        scene = Scene()
        draw_atoms_in_spacefilling_style(scene, [], Look.CARTOON, Color(255, 255, 255), 0.05)
        self.assertEqual(len(scene.nodes), 0)

    def test_draw_atoms_in_spacefilling_style(self):
        """Draw simple molecule in spacefilling style."""
        scene = Scene()
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        draw_atoms_in_spacefilling_style(scene, atoms, Look.CARTOON, Color(255, 255, 255), 0.05)
        self.assertEqual(len(scene.nodes), 4)


class TestDrawBondsInTubeStyle(unittest.TestCase):
    """Test the draw_bonds_in_tube_style function."""

    def test_draw_bonds_in_tube_style_without_supplying_bonds(self):
        """Draw simple molecule in tube style without supplying bonds."""
        scene = Scene()
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        draw_bonds_in_tube_style(
            scene,
            atoms,
            [],
            Style.TUBE,
            Look.CARTOON,
            CylinderCapType.NO_CAP,
            Color(255, 255, 255),
            0.05,
        )
        self.assertEqual(len(scene.nodes), 0)

    def test_draw_bonds_in_tube_style_without_supplying_atoms(self):
        """Draw simple molecule in tube style without supplying atoms throws KeyError."""
        scene = Scene()
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        with self.assertRaises(KeyError):
            draw_bonds_in_tube_style(
                scene,
                [],
                bonds,
                Style.TUBE,
                Look.CARTOON,
                CylinderCapType.NO_CAP,
                Color(255, 255, 255),
                0.05,
            )

    def test_draw_bonds_in_tube_style(self):
        """Draw simple molecule in tube style."""
        scene = Scene()
        atoms = [
            Atom(0, "C", (0, 0, 0)),
            Atom(1, "H", (1, 0, 0)),
            Atom(2, "H", (0, 1, 0)),
            Atom(3, "H", (0, 0, 1)),
        ]
        bonds = [
            Bond(0, 1, 1),
            Bond(0, 2, 1),
            Bond(0, 3, 1),
        ]
        draw_bonds_in_tube_style(
            scene,
            atoms,
            bonds,
            Style.TUBE,
            Look.CARTOON,
            CylinderCapType.NO_CAP,
            Color(255, 255, 255),
            0.05,
        )
        self.assertEqual(len(scene.nodes), 6)  # 3 x 2 half bonds.


class TestDrawMolecule(unittest.TestCase):
    """Test the draw_molecule function."""

    def test_draw_molecule_in_spacefilling_style(self):
        """Draw simple molecule in spacefilling style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "H", (1, 0, 0)),
                Atom(2, "H", (0, 1, 0)),
                Atom(3, "H", (0, 0, 1)),
            ],
            [
                Bond(0, 1, 1),
                Bond(0, 2, 1),
                Bond(0, 3, 1),
            ],
            Style.SPACEFILLING,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 4)
        self.assertEqual(len(svg.fills), 4)

    def test_draw_molecule_with_only_single_bonds_in_ball_and_stick_style(self):
        """Draw simple molecule in ball and stick style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "H", (1, 0, 0)),
                Atom(2, "H", (0, 1, 0)),
                Atom(3, "H", (0, 0, 1)),
            ],
            [
                Bond(0, 1, 1),
                Bond(0, 2, 1),
                Bond(0, 3, 1),
            ],
            Style.BALL_AND_STICK,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 10)
        self.assertEqual(len(svg.fills), 10)

    def test_draw_molecule_with_single_and_triple_bonds_in_ball_and_stick_style(self):
        """Draw simple molecule in ball and stick style without throwing DivideByZeroError."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "H", (1, 0, 0)),
                Atom(2, "N", (0, 1, 0)),
            ],
            [
                Bond(0, 1, 1),
                Bond(0, 2, 3),
            ],
            Style.BALL_AND_STICK,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 11)
        self.assertEqual(len(svg.fills), 11)

    def test_draw_molecule_with_multi_colored_triple_bond_in_ball_and_stick_style(self):
        """Draw simple molecule with multi-colored triple bond in ball and stick style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "N", (1, 0, 0)),
            ],
            [
                Bond(0, 1, 3),
            ],
            Style.BALL_AND_STICK,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 8)
        self.assertEqual(len(svg.fills), 8)

    def test_draw_molecule_with_single_colored_triple_bond_in_ball_and_stick_style(self):
        """Draw simple molecule with single-colored triple bond in ball and stick style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "C", (1, 0, 0)),
            ],
            [
                Bond(0, 1, 3),
            ],
            Style.BALL_AND_STICK,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 5)
        self.assertEqual(len(svg.fills), 5)

    def test_draw_molecule_in_wireframe_style(self):
        """Draw simple molecule in wireframe style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "H", (1, 0, 0)),
                Atom(2, "H", (0, 1, 0)),
                Atom(3, "H", (0, 0, 1)),
            ],
            [
                Bond(0, 1, 1),
                Bond(0, 2, 1),
                Bond(0, 3, 1),
            ],
            Style.WIREFRAME,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 6)
        self.assertEqual(len(svg.fills), 6)

    def test_draw_molecule_in_tube_style(self):
        """Draw simple molecule in tube style."""
        svg = draw_molecule(
            [
                Atom(0, "C", (0, 0, 0)),
                Atom(1, "H", (1, 0, 0)),
                Atom(2, "H", (0, 1, 0)),
                Atom(3, "H", (0, 0, 1)),
            ],
            [
                Bond(0, 1, 1),
                Bond(0, 2, 1),
                Bond(0, 3, 1),
            ],
            Style.TUBE,
            Look.CARTOON,
            50,
        )
        self.assertEqual(len(svg.objects), 6)
        self.assertEqual(len(svg.fills), 6)
