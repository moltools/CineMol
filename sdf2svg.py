#!/usr/bin/env python3
import argparse

import cinemol_py.parsing as parsing
import cinemol_py.styles as styles 
import cinemol_py.types as types 
import cinemol_py.drawing as drawing


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", required=True, type=str)
    parser.add_argument("--hs", action="store_true")
    parser.add_argument("--depiction", type=int, choices=[0, 1, 2], default=1)

    parser.add_argument("--xrot", type=float, default=0.0)
    parser.add_argument("--yrot", type=float, default=0.0)
    parser.add_argument("--zrot", type=float, default=0.0)

    return parser.parse_args()


def main() -> None:

    args = cli()
    
    view_box = None 
    options = drawing.DrawOptions_get_init()
    options.Depiction = types.Depiction(args.depiction)

    options.ShowHydrogenAtoms = args.hs 
    rotation = types.Rotation(args.xrot, args.yrot, args.zrot)
    zoom = types.Zoom(1.0)

    # atoms = empty()
    # atom = types.create_atom(0, styles.AtomType(7), types.Point3D(0.0, 0.0, 0.0), 1.0)
    # atoms = append(atoms, singleton(atom))
    # mol = types.Molecule(atoms, empty())

    with open(args.sdf, "r") as fo:
        mol_block = fo.read()
        mol = parsing.parse_sdf(mol_block)[0]
        svg, _ = drawing.draw(view_box, options, rotation, zoom, mol)
        # with open("test.svg", "w") as fo: fo.write(svg)

        print(svg)

    exit(0)


if __name__ == "__main__":
    main()
