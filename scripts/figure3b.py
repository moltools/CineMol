#!/usr/bin/env python
"""
Description:    Generate two conformers, align them and draw them.
Usage:          python figure3b.py
"""
import argparse
import typing as ty
from enum import Enum

from rdkit import Chem 
from rdkit.Chem import AllChem

from cinemol.chemistry import Style, Look, Atom, Bond, draw_molecule

def cli() -> argparse.Namespace:
    """
    Command-line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=str, required=True, help="Path to output SVG file.")
    return parser.parse_args()

def generate_conformer(mol: Chem.Mol, num_confs: int) -> Chem.Mol:
    """
    Generate a conformer of a molecule.

    :param Chem.Mol mol: Molecule to generate a conformer of.
    :param int num_confs: Number of conformers to generate.
    :return: Molecule with a generated conformer.
    :rtype: Chem.Mol
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(mol) # MMFF94
    return mol

def find_substructure(mol: Chem.Mol, smarts: str) -> ty.List[int]:
    """
    Find a substructure in a molecule.

    :param Chem.Mol mol: Molecule to find a substructure in.
    :param str smarts: SMARTS string to use for substructure search.
    :return: List of atom indices that match the substructure.
    :rtype: ty.List[int]
    """
    substructure = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(substructure)
    return [atom_index for match in matches for atom_index in match]

def main() -> None:
    """
    Driver function.
    """
    args = cli()

    # Parse molecule.
    mol = Chem.MolFromSmiles(r"C[C@@H]1[C@@H](O[C@@](C([C@H]1OC)(C)C)([C@@H](C)C[C@H](C)[C@@H]([C@@]2([C@H](O2)[C@@H](C)C=C(C)C)C)O)O)CC(=O)O")
    
    backbone_smarts = r"C1C(OC(CC1)(CCCCCCCC=C(C)C))CC"
    matched = find_substructure(mol, backbone_smarts)
    
    mol = generate_conformer(mol, 50)
    confs = [conf for conf in mol.GetConformers()]

    # Align conformers.
    AllChem.AlignMolConformers(mol, atomIds=range(mol.GetNumAtoms()), RMSlist=[0, 1])

    # Draw conformers.
    atoms, bonds = [], []

    offset = 0
    for conf in confs:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H": continue
            atom_index = atom.GetIdx()
            atoms.append(Atom(atom_index + offset, atom.GetSymbol(), conf.GetAtomPosition(atom_index)))

        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetSymbol() == "H" or bond.GetEndAtom().GetSymbol() == "H": continue
            start_index = bond.GetBeginAtomIdx()
            end_index = bond.GetEndAtomIdx()
            if start_index in matched and end_index in matched: color = (230,  25,  75)
            else: color = (120, 120, 120)
            bonds.append(Bond(start_index + offset, end_index + offset, int(bond.GetBondTypeAsDouble()), color=color, opacity=0.5))

        offset = len(atoms)

    # Draw conformers.
    svg_str = draw_molecule(atoms, bonds, Style.Wireframe, Look.Glossy, 50, verbose=True, scale=10.0)

    # Write SVG file.
    with open(args.o, "w") as f:
        f.write(svg_str)    

    exit(0)

if __name__ == "__main__":
    main()