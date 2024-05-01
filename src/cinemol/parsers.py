# -*- coding: utf-8 -*-

"""This module contains functions for parsing molecular file formats."""

import typing as ty

from cinemol.api import Atom, Bond


def parse_sdf(src: str, include_hs: bool = True) -> ty.Tuple[ty.List[Atom], ty.List[Bond]]:
    """Parse first molecule from SDF file format.

    :param src: SDF file content.
    :type src: str
    :param include_hs: Include hydrogens.
    :type include_hs: bool
    :return: Atoms and bonds.
    :rtype: ty.Tuple[ty.List[Atom], ty.List[Bond]]
    """
    atoms, bonds = [], []

    lines = src.split("\n")

    counts_line = lines[3]  # Counts line of the first molecule in the SDF file.

    atom_count = int(counts_line[0:3])
    bond_count = int(counts_line[3:6])

    atom_lines = lines[4 : 4 + atom_count]
    bond_lines = lines[4 + atom_count : 4 + atom_count + bond_count]

    atom_index = 0

    # Parse atom line.
    for atom_line in atom_lines:
        atom_index += 1

        x = float(atom_line[0:10].strip())
        y = float(atom_line[10:20].strip())
        z = float(atom_line[20:30].strip())
        atom_symbol = atom_line[31:34].strip()

        atoms.append(Atom(atom_index, atom_symbol, (x, y, z)))

    # Parse bond line.
    for bond_line in bond_lines:

        start_index = int(bond_line[0:3])
        stop_index = int(bond_line[3:6])
        bond_order = int(bond_line[6:9])

        bonds.append(Bond(int(start_index), int(stop_index), int(bond_order)))

    atom_map = {atom.index: atom for atom in atoms}
    if not include_hs:
        atoms = [atom for atom in atoms if atom.symbol != "H"]
        bonds = [
            bond
            for bond in bonds
            if (atom_map[bond.start_index].symbol != "H" and atom_map[bond.end_index].symbol != "H")
        ]

    return atoms, bonds
