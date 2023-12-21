#/usr/bin/env python3
import argparse 

from rdkit import Chem

from cinemol_py import api 

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-smi", type=str, required=True, help="SMILES string")
    parser.add_argument("-out", type=str, required=True, help="Path to output SVG file")
    return parser.parse_args()

def main() -> None:
    args = cli()
    
    mol = api.create_molecule(Chem.MolFromSmiles(args.smiles))
    svg_str = api.draw_molecule(mol)

    with open(args.out, "w") as f:
        f.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()



