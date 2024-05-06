#!/usr/bin/env python3
"""
Description:    Draw a protein as wireframe, with ligands as space-filling.
Usage:          python draw_protein_with_ligands.py -o model.svg
Dependencies:   biopython==1.83
"""
import argparse
import time

import Bio.PDB as PDB

from cinemol.geometry import Line3D, Point3D, Sphere
from cinemol.model import ModelSphere, ModelWire, Scene
from cinemol.style import Color
from cinemol.style import CoreyPaulingKoltungAtomColor as CPK
from cinemol.style import Glossy
from cinemol.style import PubChemAtomRadius as RADIUS


def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="Path to input '9lyz.pdb' PDB file. ")
    parser.add_argument("-o", type=str, help="Path to output SVG file.")
    return parser.parse_args()


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    parser = PDB.PDBParser(QUIET=True)
    protein_structure = parser.get_structure("protein", args.i)

    scene = Scene()

    # Parse protein from file and draw it as wireframe.
    # Parse ligand from file and draw it as space-filling.
    for model in protein_structure:
        for chain in model:
            previous_coordinates = None
            for residue in chain:

                # Protein wire.
                if PDB.is_aa(residue):
                    coordinates = Point3D(*residue["CA"].get_coord())

                    if previous_coordinates is None:
                        previous_coordinates = coordinates
                        continue

                    line = Line3D(previous_coordinates, coordinates)
                    color = Color(0, 0, 0)
                    scene.add_node(ModelWire(line, color, width=0.5, opacity=0.75))

                    previous_coordinates = coordinates

                # Ligand space-filling.
                if residue.get_resname() in ["AMU", "NAG", "AMU"]:
                    for atom in residue.get_atoms():
                        atom_coordinates = Point3D(*atom.get_coord())
                        atom_symbol = atom.get_name()[0]
                        atom_color = CPK().get_color(atom_symbol)
                        atom_radius = RADIUS().to_angstrom(atom_symbol)

                        sphere = Sphere(atom_coordinates, atom_radius)
                        depiction = Glossy(atom_color)
                        scene.add_node(ModelSphere(sphere, depiction))

    t0 = time.time()

    # Draw scene.
    svg = scene.draw(
        resolution=100,
        include_spheres=True,
        include_wires=True,
        calculate_sphere_sphere_intersections=True,
        scale=1000.0,
        rotation_over_y_axis=-2.0,
        rotation_over_z_axis=0.5,
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

    exit()


if __name__ == "__main__":
    main()
