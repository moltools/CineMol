#/usr/bin/env python3
import argparse 

from rdkit import Chem

from cinemol_py import api 

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-sdf", type=str, required=True, help="Path to SDF file containing a molecular conformer")
    parser.add_argument("-out", type=str, required=True, help="Path to output SVG file")
    return parser.parse_args()

def main() -> None:
    args = cli()
    
    conf = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    conf = Chem.AddHs(conf)
    mol = api.create_molecule(conf)
    options = api.DrawingOptions(
        model_style=api.ModelStyle.SpaceFilling,
        art_style=api.ArtStyle.Cartoon,
        resolution=200, 
        display_hydrogen_atoms=True
    )
    svg_str = api.draw_molecule(mol=mol, options=options)

    with open(args.out, "w") as f:
        f.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()



