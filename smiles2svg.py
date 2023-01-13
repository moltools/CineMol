#!/usr/bin/env python3
import argparse 

from rdkit import Chem 
from rdkit.Chem import AllChem 

import dash
from dash import dcc

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
import cinemol_py.encoding as encoding


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", type=str, required=True)
    return parser.parse_args()


def parse_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    mol = Chem.MolFromSmiles(smiles)

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True)
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
    options.Depiction = types.Depiction(2)
    options.ShowHydrogenAtoms = True

    app = dash.Dash()
    app.layout = dash.html.Div([
        dash.html.Div(id="image-output-container"),
        dcc.Slider(0, 90, marks=None, value=0, id="x_rotation"),
        dcc.Slider(0, 90, marks=None, value=0, id="y_rotation")
    ])

    @app.callback(
        dash.Output("image-output-container", "children"), 
        [
            dash.Input("x_rotation", "value"),
            dash.Input("y_rotation", "value")
        ]
    )
    def update_image(x_rotation: int, y_rotation: int):
        rotation = types.Rotation(x_rotation, y_rotation, 0)
        zoom = types.Zoom(1.0)
        svg, _ = drawing.draw(view_box, options, rotation, zoom, mol)
        encoded = encoding.to_base64string(svg)
        return dash.html.Img(id="image", src=f"data:image/svg+xml;base64,{encoded}")

    app.run_server(port=8050, host="127.0.0.1", debug=True)

    exit(0)


if __name__ == "__main__":
    main()
