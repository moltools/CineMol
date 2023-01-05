#!/usr/bin/env python3
import argparse

from rdkit import Chem 
from rdkit.Chem import Crippen

from cinemol_py.fable_modules.fable_library.list import (
    append, 
    singleton, 
    FSharpList, 
    tail, 
    empty, 
    of_array, 
    is_empty, 
    head
)


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


def parse_mol(mol):
    atoms = empty()
    bonds = empty()

    for a in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(a.GetIdx())

        atom = types.create_atom(
            a.GetIdx(), 
            parsing.try_cast_to_atom(a.GetSymbol()), 
            types.Point3D(pos.x, pos.y, pos.z), 
            1.0
        )
        atoms = append(atoms, singleton(atom))


    bonds_tot = len([b for b in mol.GetBonds()])
    for b in mol.GetBonds():
        index = b.GetIdx()
        start = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()
        bt = b.GetBondTypeAsDouble()
        if bt == 1.0:
            bondtype = types.BondType(0)
        elif bt == 2.0:
            bondtype = types.BondType(1)
        else:
            bondtype = types.BondType(0)
        bond1 = types.BondInfo(index, start, end, bondtype, 1.0, None, None)
        bond2 = types.BondInfo(index + bonds_tot, end, start, bondtype, 1.0, None, None)
        bonds = append(bonds, singleton(bond1))
        bonds = append(bonds, singleton(bond2))

    newMol = types.Molecule(atoms, bonds)

    return newMol


def main() -> None:

    args = cli()
    
    view_box = None 
    options = drawing.DrawOptions_get_init()
    options.Depiction = types.Depiction(args.depiction)

    options.ShowHydrogenAtoms = args.hs 
    rotation = types.Rotation(args.xrot, args.yrot, args.zrot)
    zoom = types.Zoom(1.0)

    with open(args.sdf, "r") as fo:
        mol_block = fo.read()
        rdkit_mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        pattern = Chem.MolFromSmiles("C1CC(=O)N1")
        magenta = types.Color((255, 0, 255), 1.0)
        m = list(rdkit_mol.GetSubstructMatches(pattern)[0])
        
        logps, mrs = zip(*Crippen.rdMolDescriptors._CalcCrippenContribs(rdkit_mol))
        red = (255, 0, 0)
        blue = (0, 0, 255)
        diffuse = styles.Color__Diffuse_5E38073B

        def val_to_color(val):
            if val < 0: 
                val = abs(val)
                r = val * 255
                g = 0 
                b = 0
                return types.Color((r + 100, g + 0, b + 0), 1.0)
            else: 
                val = abs(val)
                r = 0 
                g = 0
                b = val * 255
                return types.Color((r + 0, g + 0, b + 100), 1.0)

        mol = parse_mol(rdkit_mol)

        for a in mol.Atoms:
            a.Color = val_to_color(logps[a.Index])
        
        for b in mol.Bonds:
            b.StartColor = val_to_color(logps[b.Start])
            b.EndColor = val_to_color(logps[b.End]) 

        svg, _ = drawing.draw(view_box, options, rotation, zoom, mol)
        # with open("test.svg", "w") as fo: fo.write(svg)

        print(svg)


if __name__ == "__main__":
    main()
