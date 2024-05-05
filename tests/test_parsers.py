# -*- coding: utf-8 -*-

"""Contains unit tests for the cinemol.parsers module."""

import os
import unittest
from pathlib import Path

from cinemol.parsers import parse_sdf

FIXTURES = Path(__file__).parent / "fixtures"


class TestParsingSDF(unittest.TestCase):
    """Test the parsing of SDF files."""

    def test_parse_sdf_with_one_mol_below_100_atoms_with_end_line(self):
        """Test parsing of a single molecule with less than 100 atoms."""
        with open(os.path.join(FIXTURES, "test_input_001.sdf")) as file_open:
            src = file_open.read()
        atoms, bonds = parse_sdf(src)
        self.assertEqual(len(atoms), 41)
        self.assertEqual(len(bonds), 43)

    def test_parse_sdf_with_one_mol_below_100_atoms_without_end_line(self):
        """Test parsing of a single molecule with less than 100 atoms."""
        with open(os.path.join(FIXTURES, "test_input_002.sdf")) as file_open:
            src = file_open.read()
        atoms, bonds = parse_sdf(src)
        self.assertEqual(len(atoms), 41)
        self.assertEqual(len(bonds), 43)

    def test_parse_sdf_with_one_mol_below_100_atoms_excluding_hs(self):
        """Test parsing of a single molecule with less than 100 atoms."""
        with open(os.path.join(FIXTURES, "test_input_001.sdf")) as file_open:
            src = file_open.read()
        atoms, bonds = parse_sdf(src, include_hs=False)
        self.assertEqual(len(atoms), 23)
        self.assertEqual(len(bonds), 25)

    def test_parse_sdf_with_one_mol_above_100_atoms_without_end_line(self):
        """Test parsing of a single molecule with more than 100 atoms."""
        with open(os.path.join(FIXTURES, "test_input_003.sdf")) as file_open:
            src = file_open.read()
        atoms, bonds = parse_sdf(src)
        self.assertEqual(len(atoms), 100)
        self.assertEqual(len(bonds), 109)

    def test_parse_sdf_with_multiple_mols(self):
        """Test parsing of multiple molecules."""
        with open(os.path.join(FIXTURES, "test_input_004.sdf")) as file_open:
            src = file_open.read()
        atoms, bonds = parse_sdf(src)
        self.assertEqual(len(atoms), 41)
        self.assertEqual(len(bonds), 43)
