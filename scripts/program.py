#!/usr/bin/env python3
import argparse 

from rdkit import Chem

from cinemol.model import Scene, Node 
from cinemol.geometry import Point3D
from cinemol.style import CoreyPaulingKoltungAtomColor, PubChemAtomRadius, Style

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Input file path to SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Output file path to SVG file.")
    return parser.parse_args()

def main() -> None:
    args = cli()
    
    mol = Chem.MolFromMolFile(args.i, removeHs=False)
    pos = mol.GetConformer().GetPositions()

    scene = Scene()
    for i, atom in enumerate(mol.GetAtoms()):
        idx = atom.GetIdx()
        crd = Point3D(*pos[i])
        smb = atom.GetSymbol()
        rad = PubChemAtomRadius().to_angstrum(smb)
        col = CoreyPaulingKoltungAtomColor().get_color(smb)
        node = Node(idx, crd, rad, col)
        scene.add_node(node)

    print(scene)

    with open(args.o, "w") as file_open:
        file_open.write(scene.draw())

    exit(0)

if __name__ == "__main__":
    main()