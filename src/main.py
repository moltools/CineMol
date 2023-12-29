#!/usr/bin/env python3
import argparse 
from time import time

from rdkit import Chem

from cinemol.geometry import Point3D
from cinemol.shapes import Sphere, Cylinder, CapType
from cinemol.model import Scene, ModelSphere, ModelCylinder
from cinemol.style import CoreyPaulingKoltungAtomColor, PubChemAtomRadius, FillStyleType

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Output file path to SVG file.")
    parser.add_argument("-res", type=int, default=30, help="Resolution of SVG model.")
    return parser.parse_args()

def main() -> None:
    args = cli()
    
    mol = Chem.MolFromMolFile(args.i, removeHs=False)
    pos = mol.GetConformer().GetPositions()

    scene = Scene()
    fill_style = FillStyleType.Cartoon

    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        atom_color = CoreyPaulingKoltungAtomColor().get_color(atom_symbol)
        atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol) / 3

        node = ModelSphere(
            geometry=Sphere(Point3D(*pos[i]), atom_radius),
            fill_color=atom_color,
            fill_style=fill_style
        )
        scene.add_node(node) 

    for i, bond in enumerate(mol.GetBonds()):
        start_atom = bond.GetBeginAtom()
        start_index = start_atom.GetIdx()   
        start_pos = Point3D(*pos[start_index])
        start_color = CoreyPaulingKoltungAtomColor().get_color(start_atom.GetSymbol())

        end_atom = bond.GetEndAtom()
        end_index = end_atom.GetIdx()
        end_pos = Point3D(*pos[end_index])
        end_color = CoreyPaulingKoltungAtomColor().get_color(end_atom.GetSymbol())

        middle_pos = Point3D(
            (start_pos.x + end_pos.x) / 2,
            (start_pos.y + end_pos.y) / 2,
            (start_pos.z + end_pos.z) / 2
        )

        bond_radius = 0.3
        cap_style = CapType.NoCap

        start_edge = ModelCylinder(
            geometry=Cylinder(start_pos, middle_pos, bond_radius, cap_style),
            fill_color=start_color,
            fill_style=fill_style
        )
        scene.add_node(start_edge)

        end_edge = ModelCylinder(
            geometry=Cylinder(middle_pos, end_pos, bond_radius, cap_style),
            fill_color=end_color,
            fill_style=fill_style
        )
        scene.add_node(end_edge)

    t0 = time()
    svg_str = scene.draw(res=args.res, verb=True)
    print(f"Time taken to generate SVG for {scene}: {(time() - t0) * 1000:.3f} ms")

    with open(args.o, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()