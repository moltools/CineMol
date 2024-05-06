# -*- coding: utf-8 -*-

"""Command line interface for :mod:`cinemol`."""

import logging
from time import time

import click

from cinemol.api import Look, Style, draw_molecule
from cinemol.parsers import parse_sdf

__all__ = ["main"]

logger = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.argument("input", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option(
    "-s",
    "--style",
    type=click.Choice(["spacefilling", "ballandstick", "tube", "wireframe"]),
    required=False,
    default="ballandstick",
    show_default=True,
    help="Depiction style.",
)
@click.option(
    "-l",
    "--look",
    type=click.Choice(["cartoon", "glossy"]),
    required=False,
    default="cartoon",
    show_default=True,
    help="Look of the depiction.",
)
@click.option(
    "-r",
    "--resolution",
    type=int,
    required=False,
    default=30,
    show_default=True,
    help="Resolution of SVG model.",
)
@click.option(
    "-sc",
    "--scale",
    type=float,
    required=False,
    default=1.0,
    show_default=True,
    help="Scale of the model.",
)
@click.option(
    "-fl",
    "--focal-length",
    type=float,
    required=False,
    default=None,
    show_default=True,
    help="Focal length of the camera.",
)
@click.option(
    "-rx",
    "--rotation-over-x-axis",
    type=int,
    required=False,
    default=0.0,
    show_default=True,
    help="Rotation over x-axis.",
)
@click.option(
    "-ry",
    "--rotation-over-y-axis",
    type=int,
    required=False,
    default=0.0,
    show_default=True,
    help="Rotation over y-axis.",
)
@click.option(
    "-rz",
    "--rotation-over-z-axis",
    type=int,
    required=False,
    default=0.0,
    show_default=True,
    help="Rotation over z-axis.",
)
@click.option(
    "-v",
    "--verbosity",
    type=str,
    required=False,
    default="INFO",
    show_default=True,
    help="Verbosity level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
)
@click.option("-hs", "--include-hydrogens", is_flag=True, help="Include hydrogen atoms.")
def main(
    input: str,
    output: str,
    style: str,
    look: str,
    resolution: int,
    scale: float,
    focal_length: float,
    rotation_over_x_axis: int,
    rotation_over_y_axis: int,
    rotation_over_z_axis: int,
    include_hydrogens: bool,
    verbosity: bool,
) -> None:
    """Run the CineMol command line interface."""
    logging.basicConfig(level=verbosity)

    # Map CLI arguments to API.
    mapped_style = {
        "spacefilling": Style.SPACEFILLING,
        "ballandstick": Style.BALL_AND_STICK,
        "tube": Style.TUBE,
        "wireframe": Style.WIREFRAME,
    }[style]

    # Map CLI arguments to API.
    mapped_look = {"cartoon": Look.CARTOON, "glossy": Look.GLOSSY}[look]

    # Parse SDF file.
    sdf_str = open(input, "r", encoding="utf-8").read()

    atoms, bonds = parse_sdf(sdf_str)

    # Draw molecule, timed.
    t0 = time()
    svg = draw_molecule(
        atoms=atoms,
        bonds=bonds,
        style=mapped_style,
        look=mapped_look,
        resolution=resolution,
        view_box=None,
        rotation_over_x_axis=rotation_over_x_axis,
        rotation_over_y_axis=rotation_over_y_axis,
        rotation_over_z_axis=rotation_over_z_axis,
        scale=scale,
        focal_length=focal_length,
        exclude_atoms=None if include_hydrogens else ["H"],
    )
    svg_str = svg.to_svg()

    runtime = (time() - t0) * 1000  # Runtime in milliseconds.

    logger.info(f" SVG written out to: {output}")
    logger.info(f" Time taken to generate SVG: {runtime:.3f} ms")
    logger.info(f" Size of SVG: {len(svg_str) / 1000:.3f} kb")

    with open(output, "w", encoding="utf-8") as file_open:
        file_open.write(svg_str)

    exit(0)


if __name__ == "__main__":
    main()
