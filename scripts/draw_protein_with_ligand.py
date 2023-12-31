#!/usr/bin/env python3
import argparse
import os
from time import time

from Bio import PDB
from rdkit import Chem
from rdkit.Chem import rdmolops

from cinemol.model import Scene, ModelWire, ModelSphere
from cinemol.geometry import Point3D, Line3D, Sphere
from cinemol.style import Color, CoreyPaulingKoltungAtomColor, PubChemAtomRadius, Glossy

def cli() -> argparse.Namespace:
    """
    Command line interface for this script.
    
    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input PDB file.")
    parser.add_argument("-r", "--resolution", type=int, default=50, help="Resolution of SVG model.")
    parser.add_argument("-hs", "--hydrogens", action="store_true", help="Include hydrogens.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output dir.")
    return parser.parse_args()

def main() -> None:
    args = cli()

    path_log = os.path.join(args.output, "log_draw_protein_with_ligand.txt")
    log_open = open(path_log, "w")

    scene = Scene()

    # Parse protein structure from PDB file as wireframe.
    parser = PDB.PDBParser(QUIET=True)
    protein_structure = parser.get_structure("protein", args.input)
    for model in protein_structure:
        for chain in model:
            prev_residue = None # Don't connect residues between chains.

            for residue_ind, residue in enumerate(chain):
                if PDB.is_aa(residue):
                    
                    if prev_residue is not None:
                        start = Point3D(*prev_residue["CA"].get_coord())
                        end = Point3D(*residue["CA"].get_coord())
                        node = ModelWire(Line3D(start, end), Color(0, 0, 0), 0.25, 0.5)
                        scene.add_node(node)
                    
                    prev_residue = residue

    # Parse ligand structure from PDB file as space-filling model.
    # Source: https://gist.github.com/ptosco/2ee199b219b27e01052a7a1433b3bd22 
    pdb_data = open(args.input, "r").read()
    mol = Chem.RWMol(Chem.MolFromPDBBlock(pdb_data))
    bonds_to_cleave = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    for b in mol.GetBonds()
    if b.GetBeginAtom().GetPDBResidueInfo().GetIsHeteroAtom()
    ^ b.GetEndAtom().GetPDBResidueInfo().GetIsHeteroAtom()}
    [mol.RemoveBond(*b) for b in bonds_to_cleave];
    hetatm_frags = [f for f in rdmolops.GetMolFrags(mol, asMols=True)
    if f.GetNumAtoms()
    and f.GetAtomWithIdx(0).GetPDBResidueInfo().GetIsHeteroAtom()]

    # ligand = max(hetatm_frags, key=lambda x: x.GetNumAtoms()) # Get largest fragment.

    for ligand in hetatm_frags:
        ligand_coordinates = ligand.GetConformer().GetPositions()

        # Add ligand to scene.
        for atom in ligand.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_pos = Point3D(*ligand_coordinates[atom_idx])
            atom_symbol = atom.GetSymbol()
            atom_color = CoreyPaulingKoltungAtomColor().get_color(atom_symbol)
            atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol)
            sphere = Sphere(atom_pos, atom_radius)
            depiction = Glossy(atom_color, 1.0)
            node = ModelSphere(sphere, depiction)
            scene.add_node(node)

    # Draw scene.
    t0 = time()
    svg_str = scene.draw(
        resolution=args.resolution, 
        verbose=True,
        calculate_cylinder_cylinder_intersections=False,
        calculate_sphere_sphere_intersections=True,
        calculate_sphere_cylinder_intersections=False,
    )
    runtime_ms =(time() - t0) * 1000 # Runtime in milliseconds.
    print(f"runtime (ms): {runtime_ms}")
    log_open.write(f"runtime (ms): {runtime_ms}\n")
    log_open.write(f"ligand smiles: {Chem.MolToSmiles(ligand)}")
    
    with open(os.path.join(args.output, "cinemol_protein_with_ligand.svg"), "w") as file_open:
        file_open.write(svg_str)

    log_open.close()

    exit(0)

if __name__ == "__main__":
    main()