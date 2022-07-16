#!/usr/bin/env python3
import os
import sys
import inspect

from sys import argv
from rdkit import Chem
from rdkit.Chem import AllChem

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from cinnamol import Atom, Bond, Drawer, Molecule, Point, Symbol


smiles = argv[1]
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xf00d)
AllChem.MMFFOptimizeMolecule(mol)
try: _ = mol.GetConformer()
except: print("No conformer found")
conformer = mol.GetConformer()

atoms = []
bonds = []

for atom in mol.GetAtoms():
    pos = conformer.GetAtomPosition(atom.GetIdx())
    pos = Point(pos.x, pos.y, pos.z)
    symbol = Symbol[atom.GetSymbol()]
    atoms.append(Atom(symbol, pos))

mol = Molecule(atoms, bonds)

drawer = Drawer(off_screen=True)
drawer.add_molecule(mol)
pixels = drawer.save_fig(filename="./out/conformer.png", window_size=(30, 30))
print(pixels.shape)
print(0)
