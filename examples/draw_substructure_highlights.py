#!/usr/bin/env python
"""
Description:    Generate and draw a daptomycin conformer and highlight individual amino acids.
Usage:          python draw_substructure_highlights.py -o model.svg
Dependencies:   rdkit==2023.9.6
"""
import argparse
import time
import typing as ty
from enum import Enum

from rdkit import Chem
from rdkit.Chem import AllChem

from cinemol.api import Atom, Bond, Look, Style, draw_molecule


class Palette(Enum):
    """
    Palette of colors.
    """

    Red = (230, 25, 75)
    Blue = (0, 130, 200)
    Green = (60, 180, 75)
    Maroon = (128, 0, 0)
    Brown = (170, 110, 40)
    Olive = (128, 128, 0)
    Teal = (0, 128, 128)
    Navy = (0, 0, 128)
    Orange = (245, 130, 48)
    Yellow = (255, 225, 25)
    Lime = (210, 245, 60)
    Cyan = (70, 240, 240)
    Purple = (145, 30, 180)
    Magenta = (240, 50, 230)
    Pink = (255, 190, 212)
    Apricot = (255, 215, 180)
    Beige = (255, 250, 200)
    Mint = (170, 255, 195)
    Lavender = (220, 190, 255)


def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=str, help="Path to output SVG file.")
    return parser.parse_args()


def generate_conformer(mol: Chem.Mol) -> Chem.Mol:
    """
    Generate a conformer of a molecule.

    :param Chem.Mol mol: Molecule to generate a conformer of.
    :return: Molecule with a generated conformer.
    :rtype: Chem.Mol
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xF00D)
    AllChem.MMFFOptimizeMolecule(mol)  # MMFF94
    return mol


def find_substructure(mol: Chem.Mol, smarts: str) -> ty.List[int]:
    """
    Find a substructure in a molecule.

    :param Chem.Mol mol: Molecule to find a substructure in.
    :param str smarts: SMARTS string to use for substructure search.
    :return: List of atom indices that match the substructure.
    :rtype: ty.List[int]
    """
    substructure = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(substructure)
    return [atom_index for match in matches for atom_index in match]


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    # Generate a conformer for daptomycin.
    mol = Chem.MolFromSmiles(
        r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C"
        r"@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H]("
        r"NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3="
        r"O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    )
    mol = generate_conformer(mol)
    pos = mol.GetConformer().GetPositions()

    # Get atom colors based on substructure matches. Last match is used for each atom.
    atom_colors = {}
    substructures = {
        "Kyn": r"[NH2][c]1[cH][cH][cH][cH][c]1[C](=[O])[CH2][CH]([N])[C](=[O])[O]",
        "Trp": r"[N][CH]([CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12)[C](=[O])",
        "3Me-Glu": r"[CH3][CH]([CH2][C](=[O])[OH])[CH]([N])[C](=[O])",
        "Asp": r"[N][CH]([CH2][C](=[O])[OH])[C](=[O])",
        "Orn": r"[NH2][CH2][CH2][CH2][CH]([N])[C](=[O])",
        "Asn": r"[NH2][C](=[O])[CH2][CH]([N])[C](=[O])",
        "Thr": r"[CH3][CH]([O])[CH]([N])[C](=[O])",
        "Ser": r"[N][CH]([CH2][OH])[C](=[O])",
        "Ala": r"[CH3][CH]([N])[C](=[O])",
        "Gly": r"[N][CH2][C](=[O])",
    }
    palette = list(Palette)
    for i, (_, smarts) in enumerate(substructures.items()):
        for atom_index in find_substructure(mol, smarts):
            atom_colors[atom_index] = palette[i % len(palette)].value

    # Parse atoms and bonds from molecule.
    atoms, bonds = [], []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        color = atom_colors.get(atom.GetIdx(), (120, 120, 120))
        atoms.append(Atom(atom.GetIdx(), atom.GetSymbol(), pos[atom.GetIdx()], color=color))

    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol() == "H" or bond.GetEndAtom().GetSymbol() == "H":
            continue
        start_index, end_index = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds.append(Bond(start_index, end_index, int(bond.GetBondTypeAsDouble())))

    t0 = time.time()

    # Draw molecule.
    svg = draw_molecule(
        atoms=atoms,
        bonds=bonds,
        style=Style.TUBE,
        look=Look.GLOSSY,
        resolution=50,
        rotation_over_y_axis=-2.0,
        scale=1000.0,
    )

    svg_str = svg.to_svg()

    # Time in milliseconds.
    print(f"Runtime: {1000 * (time.time() - t0)} ms")

    # Get file size in kb.
    svg_size = len(svg.to_svg()) / 1000
    print(f"File size: {svg_size} kb")

    # Write SVG to file.
    with open(args.o, "w") as f:
        f.write(svg_str)


if __name__ == "__main__":
    main()
