import typing as ty 
from dataclasses import dataclass 

from rdkit import Chem 

from .fable_modules.fable_library.list import FSharpList

from .cine_mol_style import (
    AtomColorStyle__Color_Z3BD3A37D, 
    AtomColorStyle
)
from .cine_mol_types import (
    Chem_AtomType_FromString_Z721C83C5,
    Chem_BondType_FromString_Z721C83C5,
    Chem_Atom, 
    Chem_Bond, 
    Chem_Molecule,
    Color,
    Geometry_Point3D,
    Drawing_DrawingOptions_New
)
from .cine_mol_drawing import draw 
    
def create_atom(acc: int, atom: Chem.Atom) -> Chem_Atom:
    atom_symbol = atom.GetSymbol()
    pos = Geometry_Point3D(X=0.0, Y=0.0, Z=0.0)

    return Chem_Atom(
        Index=acc,
        Type=Chem_AtomType_FromString_Z721C83C5(atom_symbol),
        Color=AtomColorStyle__Color_Z3BD3A37D(AtomColorStyle(0), Chem_AtomType_FromString_Z721C83C5(atom_symbol)),
        Opacity=1.0,
        Position=pos,
        Radius=1.0
    )

def create_bond(acc: int, bond: Chem.Bond) -> Chem.Bond:
    bond_type_str = str(bond.GetBondType())
    
    return Chem_Bond(
        BeginIndex=acc,
        EndIndex=acc + 1,
        Type=Chem_BondType_FromString_Z721C83C5(bond_type_str), 
        BeginAtomIndex=bond.GetBeginAtomIdx(),
        EndAtomIndex=bond.GetEndAtomIdx(),
        Opacity=None,
        Color=None, 
        Radius=0.5
    )

def create_molecule(mol: Chem.Mol) -> Chem_Molecule:
    atoms = FSharpList(Chem_Atom, None)
    bonds = FSharpList(Chem_Bond, None)

    num_items = 0
        
    for atom in mol.GetAtoms():
        parsed = create_atom(num_items, atom)
        atoms = FSharpList(parsed, atoms)
        num_items += 1

    for bond in mol.GetBonds():
        parsed = create_bond(num_items, bond)
        bonds = FSharpList(parsed, bonds)
        num_items += 2

    return Chem_Molecule(Atoms=atoms, Bonds=bonds)

def draw_molecule(mol: Chem_Molecule) -> str:
    options = Drawing_DrawingOptions_New()
    return str(draw(mol, options))
