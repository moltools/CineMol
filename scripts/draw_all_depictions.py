#!/usr/bin/env python3
import argparse 
import os 
import subprocess
from time import time

def cli() -> argparse.Namespace:
    """
    Command line interface for this script.
    
    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input SDF file.")
    parser.add_argument("-r", "--resolution", type=int, default=50, help="Resolution of SVG model.")
    parser.add_argument("-hs", "--hydrogens", action="store_true", help="Include hydrogens.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def run_command(command: str) -> None:
    """
    Run a command in the shell.
    
    :param command: Command to run.
    :type command: str
    """
    subprocess.run(command, shell=True)

def main() -> None:
    args = cli()

    hs = "-hs" if args.hydrogens else ""
    commands = [
        # Glossy depictions.
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_glossy_spacefilling.svg')} -r {args.resolution} -d spacefilling -l glossy {hs}",
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_glossy_ballandstick.svg')} -r {args.resolution} -d ballandstick -l glossy {hs}",
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_glossy_tube.svg')} -r {args.resolution} -d tube -l glossy {hs}",

        # Cartoon depictions.
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_cartoon_spacefilling.svg')} -r {args.resolution} -d spacefilling -l cartoon {hs}",
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_cartoon_ballandstick.svg')} -r {args.resolution} -d ballandstick -l cartoon {hs}",
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_cartoon_tube.svg')} -r {args.resolution} -d tube -l cartoon {hs}",

        # Wireframe.
        f"cinemol -i {args.input} -o {os.path.join(args.output, 'cinemol_wireframe.svg')} -d wireframe {hs}",
    ]

    path_log = os.path.join(args.output, "log_draw_all_depictions.txt")
    log_open = open(path_log, "w")
    
    for i, command in enumerate(commands):
        t0 = time()
        run_command(command)
        runtime_ms =(time() - t0) * 1000 # Runtime in milliseconds.
        log_open.write(f"command: '{command}'\nruntime (ms): {runtime_ms}\n\n")
        padding = len(str(len(commands)))
        print(f"{i}".zfill(padding) + f"/{len(commands)}", end="\r")

    log_open.close()
    
    exit(0)

if __name__ == "__main__":
    main()