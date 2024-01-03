#!/usr/bin/env python3
"""
Description:       Main entry point for CineMol command line interface.
Usage:             cinemol -h
"""
import argparse 
from time import time

from cinemol.version import version
from cinemol.chemistry import Style, Look, draw_molecule
from cinemol.parsers import parse_sdf

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
    parser.add_argument("-r", type=int, default=30, help="Resolution of SVG model (default: 50).")
    parser.add_argument("-sc", type=float, default=1.0, help="Scale of the model (default: 1.0).")
    parser.add_argument("-rx", type=int, default=0.0, help="Rotation over x-axis (default: 0.0).")
    parser.add_argument("-ry", type=int, default=0.0, help="Rotation over y-axis (default: 0.0).")
    parser.add_argument("-rz", type=int, default=0.0, help="Rotation over z-axis (default: 0.0).")
    parser.add_argument("--hs", action="store_true", help="Include hydrogens (default: False).")
    parser.add_argument("--vb", action="store_true", help="Verbose mode (default: False).")
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

    # Parse SDF file.
    sdf_str = open(args.i, "r").read()
    atoms, bonds = parse_sdf(sdf_str, args.hs)
    
    # Draw molecule, timed.
    t0 = time()
    svg_str = draw_molecule(atoms, bonds, args.s, args.l, args.r, args.rx, args.ry, args.rz, args.vb, args.sc)
    runtime = (time() - t0) * 1000 # Runtime in milliseconds.

    if args.vb:
        print(f"SVG written out to: {args.o}")
        print(f"> Time taken to generate SVG: {runtime:.3f} ms")
        print(f"> Size of SVG: {len(svg_str) / 1000:.3f} kb")


    with open(args.o, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()