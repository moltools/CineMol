"""
Draw a molecule using CineMol and RDKit.
"""
from rdkit import Chem

from cinemol.geometry import Point3D, Sphere, Cylinder, CylinderCapType 
from cinemol.model import Scene, ModelSphere, ModelCylinder
from cinemol.style import Depiction, Cartoon, Glossy, CoreyPaulingKoltungAtomColor, PubChemAtomRadius

from .depiction import Style

def draw_molecule(mol: Chem.Mol, style: Style, resolution: int) -> str:
    """
    Draw a molecule using CineMol and RDKit.

    :param Chem.Mol mol: The molecule to draw.
    :param Style style: The style of the depiction.
    :param int resolution: The resolution of the SVG model.
    :return: The SVG string.
    :rtype: str
    """
    print(f"Currently argument 'style={style}' is ignored.")

    pos = mol.GetConformer().GetPositions()

    scene = Scene()

    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        atom_color = CoreyPaulingKoltungAtomColor().get_color(atom_symbol)
        atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol) / 3

        node = ModelSphere(
            geometry=Sphere(Point3D(*pos[i]), atom_radius),
            depiction=Cartoon(atom_color)
        )
        scene.add_node(node) 

    for i, bond in enumerate(mol.GetBonds()):
        start_atom = bond.GetBeginAtom()
        start_index = start_atom.GetIdx()   
        start_pos = Point3D(*pos[start_index])
        start_color = CoreyPaulingKoltungAtomColor().get_color(start_atom.GetSymbol())

        end_atom = bond.GetEndAtom()
        end_index = end_atom.GetIdx()
        end_pos = Point3D(*pos[end_index])
        end_color = CoreyPaulingKoltungAtomColor().get_color(end_atom.GetSymbol())

        middle_pos = Point3D(
            (start_pos.x + end_pos.x) / 2,
            (start_pos.y + end_pos.y) / 2,
            (start_pos.z + end_pos.z) / 2
        )

        bond_radius = 0.3
        cap_style = CylinderCapType.NoCap

        start_edge = ModelCylinder(
            geometry=Cylinder(start_pos, middle_pos, bond_radius, cap_style),
            depiction=Cartoon(start_color)
        )
        scene.add_node(start_edge)

        end_edge = ModelCylinder(
            geometry=Cylinder(middle_pos, end_pos, bond_radius, cap_style),
            depiction=Cartoon(end_color)
        )
        scene.add_node(end_edge)

    svg_str = scene.draw(resolution=resolution, verbose=True)

    return svg_str