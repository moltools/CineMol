#!/usr/bin/env python3
import argparse 
from time import time

from rdkit import Chem

from cinemol_rdkit.depiction import Style
from cinemol_rdkit.draw_molecule import draw_molecule

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draw a ball-and-stick represention of a  molecule using CineMol and RDKit.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file path to SVG file.")
    parser.add_argument("-res", "--resolution", type=int, default=30, help="Resolution of SVG model.")
    return parser.parse_args()

def main() -> None:
    args = cli()
    mol = Chem.MolFromMolFile(args.input, removeHs=False)

    t0 = time()
    svg_str = draw_molecule(mol, Style.BallAndStick, args.resolution)
    print(f"Time taken to generate SVG: {(time() - t0) * 1000:.3f} ms")

    with open(args.output, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()