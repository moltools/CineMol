#/usr/bin/env python3
import argparse 

from rdkit import Chem
from Bio import PDB 


from cinemol_py import api 

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", type=str, required=True, help="Path to PDB file containing a protein with/without ligands")
    parser.add_argument("-out", type=str, required=True, help="Path to output SVG file")
    return parser.parse_args()

def main() -> None:
    args = cli()

    parser = PDB.PDBParser(QUIET=True)
    protein_structure = parser.get_structure('protein', args.pdb)
    
    for model in protein_structure:
        for chain in model:
            for idx, residue in enumerate(chain):
                if PDB.is_aa(residue): # Check if the residue is an amino acid (standard amino acids)
                    pos = residue['CA'].get_coord()
                    radius = 0.67 * 5
                    print(residue.get_resname(), pos, radius)

    ligand = parser.get_structure('ligand', args.pdb)
    for atom in ligand.get_atoms():
        print(dir(atom))
    # ligand = Chem.AddHs(ligand)
    # atoms = [atom for atom in ligand.GetAtoms()]
    # print(len(atoms))

    # TODO: get penicillin G ligand as structure

    # options = api.DrawingOptions(
    #     model_style=api.ModelStyle.SpaceFilling,
    #     art_style=api.ArtStyle.Cartoon,
    #     resolution=200, 
    #     display_hydrogen_atoms=True
    # )
    # svg_str = api.draw_molecule(mol=mol, options=options)

    # with open(args.out, "w") as f:
    #     f.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()



