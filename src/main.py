#!/usr/bin/env python3
import argparse 
from time import time

from rdkit import Chem

from cinemol.version import version
from cinemol_rdkit.depiction import Style, Look 
from cinemol_rdkit.draw_molecule import draw_molecule

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draw a ball-and-stick represention of a  molecule using CineMol and RDKit.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file path to SVG file.")
    parser.add_argument("-d", "--depiction", type=str, required=False, default="ballandstick", choices=["spacefilling", "ballandstick", "tube", "wireframe"], help="Depiction style.")
    parser.add_argument("-l", "--look", type=str, required=False, default="cartoon", choices=["cartoon", "glossy"], help="Look of the depiction.")
    parser.add_argument("-r", "--resolution", type=int, default=50, help="Resolution of SVG model.")
    parser.add_argument("-rx", "--rotax", type=int, default=0.0, help="Rotation over x-axis.")
    parser.add_argument("-ry", "--rotay", type=int, default=0.0, help="Rotation over y-axis.")
    parser.add_argument("-rz", "--rotaz", type=int, default=0.0, help="Rotation over z-axis.")
    parser.add_argument("-hs", "--hydrogens", action="store_true", help="Include hydrogens.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode.")
    parser.add_argument("--version", action="version", version=f"CineMol {version}")
    return parser.parse_args()

def main() -> None:
    args = cli()
    mol = Chem.MolFromMolFile(args.input, removeHs=not args.hydrogens)
    Chem.Kekulize(mol)

    if args.depiction == "spacefilling":
        depiction = Style.SpaceFilling
    elif args.depiction == "ballandstick":
        depiction = Style.BallAndStick
    elif args.depiction == "tube":
        depiction = Style.Tube
    elif args.depiction == "wireframe":
        depiction = Style.Wireframe
    else:
        raise ValueError(f"Unknown depiction style: '{args.depiction}'")
    
    if args.look == "cartoon":
        look = Look.Cartoon
    elif args.look == "glossy":
        look = Look.Glossy
    else:
        raise ValueError(f"Unknown look: '{args.look}'")

    t0 = time()

    svg_str = draw_molecule(
        mol=mol, 
        style=depiction, 
        look=look,
        resolution=args.resolution, 
        rotation_over_x_axis=args.rotax,
        rotation_over_y_axis=args.rotay,
        rotation_over_z_axis=args.rotaz,
        verbose=args.verbose
    )

    if args.verbose:
        print(f"Time taken to generate SVG: {(time() - t0) * 1000:.3f} ms")

    with open(args.output, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()