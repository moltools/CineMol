#!/usr/bin/env python
"""
Description:    Generate multiple conformers, align three of them and draw.
Usage:          python draw_superimposed_conformers.py -o model.svg
Dependencies:   rdkit==2023.9.6
"""
import argparse
import time
import typing as ty

from rdkit import Chem
from rdkit.Chem import AllChem

from cinemol.geometry import (
    Cylinder,
    CylinderCapType,
    Line3D,
    Point3D,
    Sphere,
    get_perpendicular_lines,
)
from cinemol.model import ModelCylinder, ModelSphere, Scene
from cinemol.style import Color, Glossy
from cinemol.style import PubChemAtomRadius as RADIUS


def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=str, required=True, help="Path to output SVG file.")
    return parser.parse_args()


def generate_conformer(mol: Chem.Mol, num_confs: int) -> Chem.Mol:
    """
    Generate a conformer of a molecule.

    :param Chem.Mol mol: Molecule to generate a conformer of.
    :param int num_confs: Number of conformers to generate.
    :return: Molecule with a generated conformer.
    :rtype: Chem.Mol
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=0xF00D)
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

    # Parse molecule.
    mol = Chem.MolFromSmiles(r"C1=CC=C(C=C1)CC2=CC=CC=C2O")

    mol = generate_conformer(mol, 50)
    confs = [conf for conf in mol.GetConformers()]

    Chem.Kekulize(mol)

    # Align conformers.
    AllChem.AlignMolConformers(mol, atomIds=range(mol.GetNumAtoms()), RMSlist=[0, 1])

    # Draw example.
    scene = Scene()

    def draw_conformer(conf: Chem.Conformer, color: ty.Tuple[int, int, int]) -> None:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                continue

            atom_symbol = atom.GetSymbol()
            atom_coordinates = Point3D(*conf.GetAtomPosition(atom.GetIdx()))
            atom_radius = RADIUS().to_angstrom(atom_symbol) / 4.0

            sphere = Sphere(atom_coordinates, atom_radius)
            depiction = Glossy(color)
            scene.add_node(ModelSphere(sphere, depiction))

        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetSymbol() == "H" or bond.GetEndAtom().GetSymbol() == "H":
                continue

            start_index = bond.GetBeginAtomIdx()
            start_coordinates = Point3D(*conf.GetAtomPosition(start_index))
            end_index = bond.GetEndAtomIdx()
            end_coordinates = Point3D(*conf.GetAtomPosition(end_index))

            bond_order = int(bond.GetBondTypeAsDouble())
            temp_bond_radius = 0.2 / bond_order
            line = Line3D(start_coordinates, end_coordinates)
            lines = get_perpendicular_lines(line, temp_bond_radius * (bond_order + 1), bond_order)

            for line in lines:
                cylinder = Cylinder(line.start, line.end, 0.1, CylinderCapType.NO_CAP)
                depiction = Glossy(color)
                scene.add_node(ModelCylinder(cylinder, depiction))

    draw_conformer(confs[15], Color(230, 25, 75))
    draw_conformer(confs[10], Color(60, 180, 75))
    draw_conformer(confs[20], Color(0, 130, 200))

    t0 = time.time()

    # Draw example.
    svg = scene.draw(
        resolution=150,
        scale=1000.0,
        rotation_over_z_axis=-1.5,
        filter_nodes_for_intersecting=False,  # We have cyliners intersecting at other places than just start and end.
    )

    svg_str = svg.to_svg()

    # Time in milliseconds.
    print(f"Runtime: {1000 * (time.time() - t0)} ms")

    # Get file size in kb.
    svg_size = len(svg.to_svg()) / 1000
    print(f"File size: {svg_size} kb")

    # Write SVG file.
    with open(args.o, "w") as f:
        f.write(svg_str)

    exit(0)


if __name__ == "__main__":
    main()
