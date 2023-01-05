#!/usr/bin/env python3
import argparse 

from rdkit import Chem 
from rdkit.Chem import AllChem 

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
    parser.add_argument("--smiles", type=str, required=True)
    parser.add_argument("--hs", action="store_true")
    parser.add_argument("--depiction", type=int, choices=[0, 1, 2], default=1)

    parser.add_argument("--xrot", type=float, default=0.0)
    parser.add_argument("--yrot", type=float, default=0.0)
    parser.add_argument("--zrot", type=float, default=0.0)

    return parser.parse_args()


def parse_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    mol = Chem.MolFromSmiles(smiles)

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

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
    mol = parse_smiles(args.smiles)

    view_box = None 
    options = drawing.DrawOptions_get_init()
    options.Depiction = types.Depiction(args.depiction)

    options.ShowHydrogenAtoms = args.hs 
    rotation = types.Rotation(args.xrot, args.yrot, args.zrot)
    zoom = types.Zoom(1.0)

    svg, _ = drawing.draw(view_box, options, rotation, zoom, mol)
    print(svg)

    exit(0)


if __name__ == "__main__":
    main()
