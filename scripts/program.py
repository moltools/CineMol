#!/usr/bin/env python3
import argparse 
from time import time

from rdkit import Chem

from cinemol.model import Scene, Node 
from cinemol.style import CoreyPaulingKoltungAtomColor, PubChemAtomRadius, FillStyle

import numpy as np

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Output file path to SVG file.")
    return parser.parse_args()

def pick_random_style() -> FillStyle:
    options = [FillStyle.Cartoon, FillStyle.Glossy]
    return options[np.random.randint(len(options))]

def main() -> None:
    args = cli()
    
    mol = Chem.MolFromMolFile(args.i, removeHs=False)
    pos = mol.GetConformer().GetPositions()

    scene = Scene()
    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        node = Node(
            index=atom.GetIdx(), 
            center=pos[i], 
            radius=PubChemAtomRadius().to_angstrom(atom_symbol), 
            fill_color=CoreyPaulingKoltungAtomColor().get_color(atom_symbol),
            fill_style=pick_random_style()
        )
        scene.add_node(node)

    print(scene)

    t0 = time()
    svg_str = scene.draw(res=100, verb=True)
    print(f"Time taken to generate SVG: {(time() - t0) * 1000:.3f} ms")

    with open(args.o, "w") as file_open:
        file_open.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()