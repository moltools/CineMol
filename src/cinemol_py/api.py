import typing as ty 
from dataclasses import dataclass 
from enum import Enum 

from rdkit import Chem 

from .fable_modules.fable_library.list import FSharpList
from .fable_modules.fable_library.types import (Array, Union, Record, to_string)

from .cine_mol_style import (
    AtomColorStyle__Color_Z3BD3A37D, 
    AtomColorStyle,
    AtomRadius__Radius_Z3BD3A37D,
    AtomRadius
)
from .cine_mol_types import (
    Chem_AtomType_FromString_Z721C83C5,
    Chem_BondType_FromString_Z721C83C5,
    Chem_Atom, 
    Chem_Bond, 
    Chem_Molecule,
    Color,
    Geometry_Point3D,
    Drawing_DrawingOptions,
    Drawing_DrawingOptions_New,
    Svg_SVG,
    Svg_SVG__Body,
    ArtStyle as CineMolArtStyle,
    ModelStyle as CineMolModelStyle
)
from .cine_mol_drawing import draw 
    
class ModelStyle(Enum):
    SpaceFilling = 0
    BallAndStick = 1
    WireFrame = 2

class ArtStyle(Enum):   
    Cartoon = 0
    Glossy = 1

class DrawingOptions:
    def __init__(
        self, 
        model_style: ModelStyle,
        art_style: ArtStyle,
        resolution: int,
        display_hydrogen_atoms: bool
    ) -> None:
        self.model_style = model_style
        self.art_style = art_style
        self.resolution = resolution 
        self.display_hydrogen_atoms = display_hydrogen_atoms

    def _get_options(self) -> Drawing_DrawingOptions:
        options = Drawing_DrawingOptions_New()
        options.ArtStyle = CineMolArtStyle(self.art_style.value)
        options.ModelStyle = CineMolModelStyle(self.model_style.value)
        options.Resolution = self.resolution
        options.DisplayHydrogenAtoms = self.display_hydrogen_atoms
        return options
    
def create_atom(acc: int, atom: Chem.Atom, pos: Geometry_Point3D) -> Chem_Atom:
    atom_symbol = atom.GetSymbol()
    atom_type = Chem_AtomType_FromString_Z721C83C5(atom_symbol)
    atom_radius = AtomRadius__Radius_Z3BD3A37D(AtomRadius(0), atom_type)

    return Chem_Atom(
        Index=acc,
        Type=atom_type,
        Color=AtomColorStyle__Color_Z3BD3A37D(AtomColorStyle(0), atom_type),
        Opacity=1.0,
        Position=pos,
        Radius=atom_radius
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
    conf = mol.GetConformer()
    pos = conf.GetPositions()

    atoms = FSharpList(Chem_Atom, None)
    bonds = FSharpList(Chem_Bond, None)

    num_items = 0
        
    for i, atom in enumerate(mol.GetAtoms()):
        atom_pos = Geometry_Point3D(X=pos[i][0], Y=pos[i][1], Z=pos[i][2])
        parsed = create_atom(num_items, atom, atom_pos)
        atoms = FSharpList(parsed, atoms)
        num_items += 1

    for bond in mol.GetBonds():
        parsed = create_bond(num_items, bond)
        bonds = FSharpList(parsed, bonds)
        num_items += 2

    return Chem_Molecule(Atoms=atoms, Bonds=bonds)

def draw_molecule(mol: Chem_Molecule, options: DrawingOptions) -> str:
    options = options._get_options()    
    svg = draw(mol, options)[0]
    return str(svg)
