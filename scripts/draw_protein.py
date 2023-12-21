#/usr/bin/env python3
import argparse 
from enum import Enum

from rdkit import Chem
from rdkit.Chem import rdmolops
from Bio import PDB 
from copy import deepcopy

from cinemol_py import api 

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", type=str, required=True, help="Path to PDB file containing a protein with/without ligands")
    parser.add_argument("-out", type=str, required=True, help="Path to output SVG file")
    return parser.parse_args()

class AminoAcidCharge(Enum):
    ARG = 1
    LYS = 1 # pos charge
    HIS = 1
    ASP = -1 # neg charge
    GLU = -1
    CYS = -2 # non-polar
    TYR = 2 # polar
    SER = 2
    THR = 2
    ASN = 2
    GLN = 2
    GLY = -2
    PRO = -2
    ALA = -2
    ILE = -2
    LEU = -2
    MET = -2 
    PHE = -2
    TRP = -2 
    VAL = -2

def charge_to_color(charge: float) -> list:
    if charge == -2:
        # yellow
        return [255, 255, 0]
    elif charge == 2:
        # lighb blue
        return [0, 255, 255]
    elif charge == -1:
        # purple 
        return [255, 0, 255]
    elif charge == 1:
        # green
        return [0, 255, 0]
    else:
        return [120, 120, 120]

def main() -> None:
    args = cli()

    parser = PDB.PDBParser(QUIET=True)
    protein_structure = parser.get_structure('protein', args.pdb)
    
    wire_frame_colors = [ 
        [139, 0, 0],
        [0, 128, 0],
    ]

    nodes = []
    edges = []
    num_objs = 0
    for model in protein_structure:
        for chain_idx, chain in enumerate(model):
            prev_node = None
            for idx, residue in enumerate(chain):
                if PDB.is_aa(residue): # Check if the residue is an amino acid (standard amino acids)
                    pos = residue['CA'].get_coord()
                    # radius = 0.67 * 5
                    radius= 1 * 3.5

                    color = charge_to_color(AminoAcidCharge[residue.resname].value)

                    node = api.create_node(
                        num_objs, 
                        pos,
                        # wire_frame_colors[chain_idx],
                        # [120, 120, 120],
                        # [36,158,123],
                        color,
                        radius,
                        opacity=1.00
                    )
                    nodes.append(node)
                    num_objs += 1

                    if prev_node is not None:
                        edges.append(api.create_edge(
                            num_objs+1,
                            num_objs+2,
                            prev_node.Index,
                            node.Index
                        ))
                        num_objs += 2

                    prev_node = node 


    # https://gist.github.com/ptosco/2ee199b219b27e01052a7a1433b3bd22 
    pdb_data = open(args.pdb, "r").read()
    mol = Chem.RWMol(Chem.MolFromPDBBlock(pdb_data))
    bonds_to_cleave = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    for b in mol.GetBonds()
    if b.GetBeginAtom().GetPDBResidueInfo().GetIsHeteroAtom()
    ^ b.GetEndAtom().GetPDBResidueInfo().GetIsHeteroAtom()}
    [mol.RemoveBond(*b) for b in bonds_to_cleave];
    hetatm_frags = [f for f in rdmolops.GetMolFrags(mol, asMols=True)
    if f.GetNumAtoms()
    and f.GetAtomWithIdx(0).GetPDBResidueInfo().GetIsHeteroAtom()]
    # Get largest fragment
    hetatm_frag = max(hetatm_frags, key=lambda x: x.GetNumAtoms())

    num_nodes_in_protein = deepcopy(num_objs)
    pos = hetatm_frag.GetConformer().GetPositions()
    for i, atom in enumerate(hetatm_frag.GetAtoms()):

        node_pos = api.make_point(pos[i])
        node = api.create_atom(num_nodes_in_protein + atom.GetIdx(), atom, node_pos)
        nodes.append(node)
        num_objs += 1

    for bond in hetatm_frag.GetBonds():
        edge = api.create_bond(num_objs, bond, num_nodes_in_protein)
        edges.append(edge)
        num_objs += 2
        

    # mol = api.create_molecule(hetatm_frag)

    print(len(nodes), len(edges))
    obj = api.create_object(nodes, edges)

    options = api.DrawingOptions(
        model_style=api.ModelStyle.SpaceFilling,
        art_style=api.ArtStyle.Glossy,
        resolution=50, 
        display_hydrogen_atoms=False
    )
    svg_str = api.draw_molecule(mol=obj, options=options)

    with open(args.out, "w") as f:
        f.write(svg_str)

    exit(0)

if __name__ == "__main__":
    main()



