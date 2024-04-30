# -*- coding: utf-8 -*-

"""Command line interface for :mod:`cinemol`."""

import argparse
import logging
from time import time

import click

from cinemol.chemistry import Look, Style, draw_molecule
from cinemol.parsers import parse_sdf
from cinemol.version import VERSION

__all__ = ["main"]

logger = logging.getLogger(__name__)


def cli() -> argparse.Namespace:
    """Parse command line arguments.

    :return: The parsed arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Draw a ball-and-stick represention of a  molecule using CineMol and RDKit."
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input file path to SDF file. Only the first molecule in the SDF file is drawn.",
    )

    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output file path to SVG file."
    )

    parser.add_argument(
        "-s",
        "--style",
        type=str,
        required=False,
        default="ballandstick",
        choices=["spacefilling", "ballandstick", "tube", "wireframe"],
        help="Depiction style (default: ballandstick).",
    )

    parser.add_argument(
        "-l",
        "--look",
        type=str,
        required=False,
        default="cartoon",
        choices=["cartoon", "glossy"],
        help="Look of the depiction (default: cartoon).",
    )

    parser.add_argument(
        "-r", "--resolution", type=int, default=30, help="Resolution of SVG model (default: 50)."
    )

    parser.add_argument(
        "-sc", "--scale", type=float, default=1.0, help="Scale of the model (default: 1.0)."
    )

    parser.add_argument(
        "-fl",
        "--focal-length",
        type=float,
        default=10.0,
        help="Focal length of the camera (default: 10.0).",
    )

    parser.add_argument(
        "-rx",
        "--rotation-over-x-axis",
        type=int,
        default=0.0,
        help="Rotation over x-axis (default: 0.0).",
    )

    parser.add_argument(
        "-ry",
        "--rotation-over-y-axis",
        type=int,
        default=0.0,
        help="Rotation over y-axis (default: 0.0).",
    )

    parser.add_argument(
        "-rz",
        "--rotation-over-z-axis",
        type=int,
        default=0.0,
        help="Rotation over z-axis (default: 0.0).",
    )

    parser.add_argument(
        "-hs",
        "--include-hydrogens",
        action="store_true",
        help="Include hydrogens (default: False).",
    )

    parser.add_argument(
        "-vb", "--verbose", action="store_true", help="Verbose mode (default: False)."
    )

    parser.add_argument("-v", "--version", action="version", version=f"CineMol {VERSION}")

    args = parser.parse_args()

    args.style = {
        "spacefilling": Style.SPACEFILLING,
        "ballandstick": Style.BALL_AND_STICK,
        "tube": Style.TUBE,
        "wireframe": Style.WIREFRAME,
    }[args.style]

    args.look = {"cartoon": Look.CARTOON, "glossy": Look.GLOSSY}[args.look]

    return args


@click.group()
@click.version_option()
def main() -> None:
    """Run the CineMol command line interface."""
    args = cli()

    # Parse SDF file.
    sdf_str = open(args.i, "r", encoding="utf-8").read()

    atoms, bonds = parse_sdf(sdf_str)

    # Draw molecule, timed.
    t0 = time()
    svg = draw_molecule(
        atoms=atoms,
        bonds=bonds,
        style=args.style,
        look=args.look,
        resolution=args.resolution,
        view_box=None,
        rotation_over_x_axis=args.rotation_over_x_axis,
        rotation_over_y_axis=args.rotation_over_y_axis,
        rotation_over_z_axis=args.rotation_over_z_axis,
        scale=args.scale,
        focal_length=args.focal_length,
        exclude_atoms=None if args.include_hydrogens else ["H"],
        verbose=args.verbose,
    )
    svg_str = svg.to_svg()

    runtime = (time() - t0) * 1000  # Runtime in milliseconds.

    if args.vb:
        logger.info(f"SVG written out to: {args.o}")
        logger.info(f"> Time taken to generate SVG: {runtime:.3f} ms")
        logger.info(f"> Size of SVG: {len(svg_str) / 1000:.3f} kb")

    with open(args.o, "w", encoding="utf-8") as file_open:
        file_open.write(svg_str)

    exit(0)


if __name__ == "__main__":
    main()
