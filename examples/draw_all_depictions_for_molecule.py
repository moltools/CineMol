#!/usr/bin/env python3
"""
Description:    Draw space-filling, ball-and-stick, and tube representations of a molecule.
Usage:          python draw_all_depictions_for_molecule.py -i /path/to/sdf/file -o /path/to/output/dir
"""

import argparse
import os
import subprocess


def cli() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: The parsed arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="Draw space-filling, ball-and-stick, and tube representations of a molecule."
    )
    parser.add_argument("-i", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Output file path to output directory.")
    return parser.parse_args()


def run_command(command: str) -> None:
    """
    Run a command in the shell.

    :param command: The command to run.
    :type command: str
    """
    subprocess.run(command, shell=True)


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    commands = [
        "cinemol -s spacefilling -l cartoon -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'cartoon_spacefilling.svg')}",
        "cinemol -s ballandstick -l cartoon -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'cartoon_ballandstick.svg')}",
        "cinemol -s tube -l cartoon -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'cartoon_tube.svg')}",
        "cinemol -s spacefilling -l glossy -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'glossy_spacefilling.svg')}",
        "cinemol -s ballandstick -l glossy -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'glossy_ballandstick.svg')}",
        "cinemol -s tube -l glossy -r 100 -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'glossy_tube.svg')}",
        "cinemol -s wireframe -sc 1000.0 -hs"
        + f" {args.i} {os.path.join(args.o, 'wireframe.svg')}",
    ]

    for command in commands:
        print(command)
        run_command(command)

    exit(0)


if __name__ == "__main__":
    main()
