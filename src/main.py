#!/usr/bin/env python3
"""
Description:       Main entry point for CineMol command line interface.
Usage:             cinemol -h
"""
import argparse 
from time import time

from rdkit import Chem

from cinemol.version import version
from cinemol.chemistry import Style, Look, Atom, Bond, draw_molecule

def cli() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    :return: The parsed arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Draw a ball-and-stick represention of a  molecule using CineMol and RDKit.")
    parser.add_argument("-i", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Output file path to SVG file.")
    parser.add_argument("-s", type=str, required=False, default="ballandstick", choices=["spacefilling", "ballandstick", "tube", "wireframe"], help="Depiction style (default: ballandstick).")
    parser.add_argument("-l", type=str, required=False, default="cartoon", choices=["cartoon", "glossy"], help="Look of the depiction (default: cartoon).")
    parser.add_argument("-r", type=int, default=50, help="Resolution of SVG model (default: 50).")
    parser.add_argument("-rx", type=int, default=0.0, help="Rotation over x-axis (default: 0.0).")
    parser.add_argument("-ry", type=int, default=0.0, help="Rotation over y-axis (default: 0.0).")
    parser.add_argument("-rz", type=int, default=0.0, help="Rotation over z-axis (default: 0.0).")
    parser.add_argument("-hs", action="store_true", help="Include hydrogens (default: False).")
    parser.add_argument("-vb", action="store_true", help="Verbose mode (default: False).")
    parser.add_argument("-v", action="version", version=f"CineMol {version}")
    args = parser.parse_args()

    args.s = {"spacefilling": Style.SpaceFilling, "ballandstick": Style.BallAndStick, "tube": Style.Tube, "wireframe": Style.Wireframe}[args.s]
    args.l = {"cartoon": Look.Cartoon, "glossy": Look.Glossy}[args.l]

    return args 

def main() -> None:
    """
    Main entry point for CineMol command line interface.
    """
    args = cli()

    # Parse molecule from SDF file, and kekulize it.
    mol = Chem.MolFromMolFile(args.i, removeHs=not args.hs)
    pos = mol.GetConformer().GetPositions()
    Chem.Kekulize(mol)
    
    # Parse atoms.
    atoms = []
    for atom in mol.GetAtoms():
        x, y, z = pos[atom.GetIdx()]
        atoms.append(Atom(atom.GetIdx(), atom.GetSymbol(), (x, y, z)))

    # Parse bonds.
    bonds = []
    for bond in mol.GetBonds():
        bond_order = int(bond.GetBondTypeAsDouble())
        bonds.append(Bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_order))
    
    # Draw molecule, timed.
    t0 = time()
    svg_str = draw_molecule(atoms, bonds, args.s, args.l, args.r, args.rx, args.ry, args.rz, args.vb)
    runtime = (time() - t0) * 1000 # Runtime in milliseconds.

    if args.vb:
        print(f"Time taken to generate SVG: {runtime:.3f} ms")

    with open(args.o, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()