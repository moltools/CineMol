"""
Draw a molecule using CineMol and RDKit.
"""
from rdkit import Chem

from cinemol.geometry import Line3D, Point3D, Sphere, Cylinder, CylinderCapType, get_perpendicular_lines
from cinemol.model import Scene, ModelSphere, ModelCylinder, ModelWire
from cinemol.style import Color, Cartoon, Glossy, CoreyPaulingKoltungAtomColor, PubChemAtomRadius

from .depiction import Style, Look

def draw_molecule(
    mol: Chem.Mol, 
    style: Style, 
    look: Look,
    resolution: int,
    rotation_over_x_axis: float = 0.0,
    rotation_over_y_axis: float = 0.0,
    rotation_over_z_axis: float = 0.0,
    verbose: bool = False
) -> str:
    """
    Draw a molecule using CineMol and RDKit.

    :param Chem.Mol mol: The molecule to draw.
    :param Style style: The style of the depiction.
    :param int resolution: The resolution of the SVG model.
    :param float rotation_over_x_axis: The rotation over the x-axis.
    :param float rotation_over_y_axis: The rotation over the y-axis.
    :param float rotation_over_z_axis: The rotation over the z-axis.
    :param bool verbose: Whether to print verbose output.
    :return: The SVG string.
    :rtype: str
    """
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"mol must be of type Chem.Mol, not {type(mol)}")
    
    if not isinstance(style, Style):
        raise TypeError(f"style must be of type Style, not {type(style)}")
    
    if not isinstance(look, Look):
        raise TypeError(f"look must be of type Look, not {type(look)}")
    
    pos = mol.GetConformer().GetPositions()
    scene = Scene()

    # Wireframe has a separate implementation that is faster than for geometric shapes.  
    if style == Style.Wireframe:
        
        for bond in mol.GetBonds():
            start_atom = bond.GetBeginAtom()
            start_symbol = start_atom.GetSymbol()
            start_color = CoreyPaulingKoltungAtomColor().get_color(start_symbol)
            start_pos = Point3D(*pos[start_atom.GetIdx()])

            end_atom = bond.GetEndAtom()
            end_symbol = end_atom.GetSymbol()
            end_color = CoreyPaulingKoltungAtomColor().get_color(end_symbol)
            end_pos = Point3D(*pos[end_atom.GetIdx()])

            if start_color != end_color:
                middle_pos = Point3D(
                    (start_pos.x + end_pos.x) / 2,
                    (start_pos.y + end_pos.y) / 2,
                    (start_pos.z + end_pos.z) / 2
                )
                scene.add_node(ModelWire(Line3D(start_pos, middle_pos), start_color, 0.05, 1.0))
                scene.add_node(ModelWire(Line3D(middle_pos, end_pos), end_color, 0.05, 1.0))
            
            else:
                scene.add_node(ModelWire(Line3D(start_pos, end_pos), start_color, 0.05, 1.0))

        svg_str = scene.draw(
            resolution=resolution, 
            verbose=verbose, 
            rotation_over_x_axis=rotation_over_x_axis,
            rotation_over_y_axis=rotation_over_y_axis,
            rotation_over_z_axis=rotation_over_z_axis,
            include_spheres=False,
            include_cylinders=False,
            include_wires=True,
            calculate_sphere_sphere_intersections=False,
            calculate_sphere_cylinder_intersections=False,
            calculate_cylinder_sphere_intersections=False,
            calculate_cylinder_cylinder_intersections=False
        )
        return svg_str
    
    # Set default settings for depiction when not wireframe.
    include_spheres = False 
    include_cylinders = False 
    calculate_sphere_sphere_intersections = False
    calculate_sphere_cylinder_intersections = False
    calculate_cylinder_sphere_intersections = False
    calculate_cylinder_cylinder_intersections = False
    stroke_color = Color(0, 0, 0)
    stroke_width = 0.05
    opacity = 1.0
    bond_radius = 0.2
    radius_modifier = 1.0
    cap_style = CylinderCapType.NoCap 

    # Set settings for depiction.
    if style == Style.SpaceFilling:
        include_spheres = True 
        calculate_sphere_sphere_intersections = True
    
    elif style == Style.BallAndStick:
        radius_modifier = 0.25
        include_cylinders = True
        include_spheres = True 
        calculate_cylinder_sphere_intersections = True
    
    elif style == Style.Tube:
        include_cylinders = True
        cap_style = CylinderCapType.Round
        calculate_cylinder_cylinder_intersections = True

    for i, atom in enumerate(mol.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        atom_color = CoreyPaulingKoltungAtomColor().get_color(atom_symbol)
        atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol) * radius_modifier

        if look == Look.Cartoon: depiction = Cartoon(atom_color, stroke_color, stroke_width, opacity)
        elif look == Look.Glossy: depiction = Glossy(atom_color, opacity)
        else: raise ValueError(f"Unknown look: '{look}'")
        node = ModelSphere(Sphere(Point3D(*pos[i]), atom_radius), depiction)
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

        bond_type = bond.GetBondType()

        # Determine number of lines to draw for each bond.
        if style == Style.Tube:
            num_lines = 1
            temp_bond_radius = bond_radius
            lines = [Line3D(start_pos, end_pos)]

        else:
            num_lines = {
                Chem.rdchem.BondType.SINGLE: 1, 
                Chem.rdchem.BondType.DOUBLE: 2, 
                Chem.rdchem.BondType.TRIPLE: 3
            }.get(bond_type, 1)
            temp_bond_radius = bond_radius / num_lines
            lines = get_perpendicular_lines(Line3D(start_pos, end_pos), temp_bond_radius * (num_lines + 1), num_lines)

        for line in lines:
            start_pos, end_pos = line.start, line.end

            if (start_color == end_color):
                if look == Look.Cartoon: depiction = Cartoon(start_color, stroke_color, stroke_width, opacity)
                elif look == Look.Glossy: depiction = Glossy(start_color, opacity)
                edge = ModelCylinder(Cylinder(start_pos, end_pos, temp_bond_radius, cap_style), depiction)
                scene.add_node(edge)
            
            else:
                middle_pos = Point3D(
                    (start_pos.x + end_pos.x) / 2,
                    (start_pos.y + end_pos.y) / 2,
                    (start_pos.z + end_pos.z) / 2
                )

                if look == Look.Cartoon: depiction = Cartoon(start_color, stroke_color, stroke_width, opacity)
                elif look == Look.Glossy: depiction = Glossy(start_color, opacity)
                start_edge = ModelCylinder(Cylinder(start_pos, middle_pos, temp_bond_radius, cap_style), depiction)
                scene.add_node(start_edge)

                if look == Look.Cartoon: depiction = Cartoon(end_color, stroke_color, stroke_width, opacity)
                elif look == Look.Glossy: depiction = Glossy(end_color, opacity)
                end_edge = ModelCylinder(Cylinder(middle_pos, end_pos, temp_bond_radius, cap_style), depiction)
                scene.add_node(end_edge)

    svg_str = scene.draw(
        resolution=resolution, 
        verbose=verbose, 
        rotation_over_x_axis=rotation_over_x_axis,
        rotation_over_y_axis=rotation_over_y_axis,
        rotation_over_z_axis=rotation_over_z_axis,
        include_spheres=include_spheres,
        include_cylinders=include_cylinders,
        calculate_sphere_sphere_intersections=calculate_sphere_sphere_intersections,
        calculate_sphere_cylinder_intersections=calculate_sphere_cylinder_intersections,
        calculate_cylinder_sphere_intersections=calculate_cylinder_sphere_intersections,
        calculate_cylinder_cylinder_intersections=calculate_cylinder_cylinder_intersections
    )
    return svg_str