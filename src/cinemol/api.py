# -*- coding: utf-8 -*-

"""This module provides functions for drawing molecules using atoms and bonds."""

import typing as ty
from dataclasses import dataclass
from enum import Enum, auto

from cinemol.geometry import (
    Cylinder,
    CylinderCapType,
    Line3D,
    Point3D,
    Sphere,
    get_perpendicular_lines,
)
from cinemol.model import ModelCylinder, ModelSphere, ModelWire, Scene
from cinemol.style import (
    Cartoon,
    Color,
    CoreyPaulingKoltungAtomColor,
    Depiction,
    Glossy,
    PubChemAtomRadius,
)
from cinemol.svg import Svg, ViewBox


class Style(Enum):
    """The style of the depiction.

    :cvar SpaceFilling: The space-filling style.
    :cvar BallAndStick: The ball-and-stick style.
    :cvar Tube: The tube style.
    :cvar Wireframe: The wireframe style.
    """

    SPACEFILLING = auto()
    BALL_AND_STICK = auto()
    TUBE = auto()
    WIREFRAME = auto()


class Look(Enum):
    """The look of the depiction.

    :cvar Cartoon: The cartoon look.
    :cvar Glossy: The glossy look.
    """

    CARTOON = auto()
    GLOSSY = auto()


@dataclass
class Atom:
    """Represents an atom in a molecule.

    :ivar index: The index of the atom.
    :type index: int
    :ivar symbol: The symbol of the atom.
    :type symbol: str
    :ivar coordinates: The coordinates of the atom.
    :type coordinates: ty.Tuple[float, float, float]
    :ivar radius: The radius of the atom.
    :type radius: ty.Optional[float]
    :ivar color: The color of the atom.
    :type color: ty.Optional[ty.Tuple[int, int, int]]
    :ivar opacity: The opacity of the atom.
    :type opacity: float
    """

    index: int
    symbol: str
    coordinates: ty.Tuple[float, float, float]
    radius: ty.Optional[float] = None
    color: ty.Optional[ty.Tuple[int, int, int]] = None
    opacity: float = 1.0


@dataclass
class Bond:
    """A bond between two atoms in a molecule.

    :ivar start_index: The index of the start atom.
    :type start_index: int
    :ivar end_index: The index of the end atom.
    :type end_index: int
    :ivar order: The order of the bond.
    """

    start_index: int
    end_index: int
    order: int
    radius: ty.Optional[float] = None
    color: ty.Optional[ty.Tuple[int, int, int]] = None
    opacity: float = 1.0


def filter_atoms_and_bonds(
    atoms: ty.List[Atom],
    bonds: ty.List[Bond],
    exclude_atoms: ty.Optional[ty.List[str]] = None,
) -> ty.Tuple[ty.List[Atom], ty.List[Bond]]:
    """
    Filter atoms and bonds based on the given exclude atoms.

    :param atoms: The atoms in the molecule.
    :type atoms: ty.List[Atom]
    :param bonds: The bonds in the molecule.
    :type bonds: ty.List[Bond]
    :param exclude_atoms: The atoms to exclude from the depiction (list of atom symbols).
    :type exclude_atoms: ty.Optional[ty.List[str]]
    :return: The filtered atoms and bonds.
    :rtype: ty.Tuple[ty.List[Atom], ty.List[Bond]]
    """
    if exclude_atoms is not None:
        filtered_atoms, filtered_bonds, exclude_inds = [], [], []

        for atom in atoms:
            if atom.symbol not in exclude_atoms:
                filtered_atoms.append(atom)
            else:
                exclude_inds.append(atom.index)

        for bond in bonds:
            if bond.start_index not in exclude_inds and bond.end_index not in exclude_inds:
                filtered_bonds.append(bond)

        return filtered_atoms, filtered_bonds

    return atoms, bonds


def draw_bonds_in_wireframe_style(
    scene: Scene,
    atoms: ty.List[Atom],
    bonds: ty.List[Bond],
    wire_width: float = 0.05,
) -> None:
    """Draw a molecule using the given atoms and bonds in wireframe style.

    :param scene: The scene to draw the molecule in.
    :type scene: Scene
    :param atoms: The atoms in the molecule. We need this to get the position and
        color of the atoms to draw the bonds.
    :type atoms: ty.List[Atom]
    :param bonds: The bonds in the molecule.
    :type bonds: ty.List[Bond]
    :param wire_width: The width of the wires.
    :type wire_width: float
    """
    atom_map = {atom.index: atom for atom in atoms}

    for bond in bonds:
        # Get start atom specifications.
        start_atom = atom_map[bond.start_index]
        start_symbol = start_atom.symbol
        start_color = CoreyPaulingKoltungAtomColor().get_color(start_symbol)
        start_pos = Point3D(*start_atom.coordinates)

        # Get end atom specifications.
        end_atom = atom_map[bond.end_index]
        end_symbol = end_atom.symbol
        end_color = CoreyPaulingKoltungAtomColor().get_color(end_symbol)
        end_pos = Point3D(*end_atom.coordinates)

        # Determine color of bond.
        if bond.color is not None:
            start_color = Color(*bond.color)
            end_color = Color(*bond.color)
        else:
            if start_atom.color is not None:
                start_color = Color(*start_atom.color)

            if end_atom.color is not None:
                end_color = Color(*end_atom.color)

        # If colors are not the same we split the bond down the middle, and
        # draw two separate wires to represent the bond.
        if start_color != end_color:
            middle_pos = Point3D(
                (start_pos.x + end_pos.x) / 2,
                (start_pos.y + end_pos.y) / 2,
                (start_pos.z + end_pos.z) / 2,
            )
            scene.add_node(
                ModelWire(Line3D(start_pos, middle_pos), start_color, wire_width, bond.opacity)
            )
            scene.add_node(
                ModelWire(Line3D(middle_pos, end_pos), end_color, wire_width, bond.opacity)
            )

        else:
            scene.add_node(
                ModelWire(Line3D(start_pos, end_pos), start_color, wire_width, bond.opacity)
            )


def draw_atoms_in_spacefilling_style(
    scene: Scene,
    atoms: ty.List[Atom],
    look: Look,
    stroke_color: Color,
    stroke_width: float,
    radius_scale: float = 1.0,
) -> None:
    """Draw a molecule using the given atoms in space-filling style.

    :param scene: The scene to draw the molecule in.
    :type scene: Scene
    :param atoms: The atoms in the molecule.
    :type atoms: ty.List[Atom]
    :param look: The look of the depiction.
    :type look: Look
    :param stroke_color: The stroke color of the depiction.
    :type stroke_color: Color
    :param stroke_width: The stroke width of the depiction.
    :type stroke_width: float
    :param radius_scale: Divide the radius of the atoms by this value.
    :type radius_scale: float
    :raises ValueError: If an unknown look is given.
    """
    for atom in atoms:
        # Get atom specifications.
        atom_symbol = atom.symbol
        atom_color = (
            CoreyPaulingKoltungAtomColor().get_color(atom_symbol)
            if atom.color is None
            else Color(*atom.color)
        )
        atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol) * radius_scale
        atom_pos = Point3D(*atom.coordinates)

        # Determine atom look.
        if look == Look.CARTOON:
            depiction: Depiction = Cartoon(atom_color, stroke_color, stroke_width, atom.opacity)

        elif look == Look.GLOSSY:
            depiction = Glossy(atom_color, atom.opacity)

        else:
            raise ValueError(f"Unknown look: '{look}'")

        # Add atom to scene.
        scene.add_node(ModelSphere(Sphere(atom_pos, atom_radius), depiction))


def draw_bonds_in_tube_style(
    scene: Scene,
    atoms: ty.List[Atom],
    bonds: ty.List[Bond],
    tube_bond_style: Style,
    look: Look,
    cap_type: CylinderCapType,
    stroke_color: Color,
    stroke_width: float,
) -> None:
    """Draw a molecule using the given atoms and bonds in tube style.

    :param scene: The scene to draw the molecule in.
    :type scene: Scene
    :param atoms: The atoms in the molecule. We need this to get the position and
        color of the atoms to draw the bonds.
    :type atoms: ty.List[Atom]
    :param bonds: The bonds in the molecule.
    :type bonds: ty.List[Bond]
    :param tube_bond_style: The style of the depiction. Ball-and-stick will draw
        bond order as cylinders, while tube will draw any bond order as a single
        cylinder.
    :type tube_bond_style: Style
    :param look: The look of the depiction.
    :type look: Look
    :param cap_type: The cap type of the cylinders.
    :type cap_type: CylinderCapType
    :param stroke_color: The stroke color of the depiction.
    :type stroke_color: Color
    :param stroke_width: The stroke width of the depiction.
    :type stroke_width: float
    :raises ValueError: If an unknown look is given.
    """
    atom_map = {atom.index: atom for atom in atoms}

    for bond in bonds:
        # Get start atom specifications.
        start_atom = atom_map[bond.start_index]
        start_symbol = start_atom.symbol
        start_color = CoreyPaulingKoltungAtomColor().get_color(start_symbol)
        start_pos = Point3D(*start_atom.coordinates)

        # Get end atom specifications.
        end_atom = atom_map[bond.end_index]
        end_symbol = end_atom.symbol
        end_color = CoreyPaulingKoltungAtomColor().get_color(end_symbol)
        end_pos = Point3D(*end_atom.coordinates)

        # Determine color of bond.
        if bond.color is not None:
            start_color = Color(*bond.color)
            end_color = Color(*bond.color)
        else:
            if start_atom.color is not None:
                start_color = Color(*start_atom.color)

            if end_atom.color is not None:
                end_color = Color(*end_atom.color)

        # Determine number of cylinders to draw for each bond.
        # Bond order is only used for ball-and-stick style.
        bond_order = bond.order if tube_bond_style == Style.BALL_AND_STICK else 1
        bond_radius = bond.radius if bond.radius is not None else 0.2
        temp_bond_radius = bond_radius / bond_order
        line = Line3D(start_pos, end_pos)
        lines = get_perpendicular_lines(line, temp_bond_radius * (bond_order + 1), bond_order)

        # Add cylinders to scene.
        for line in lines:
            if start_color != end_color:
                middel_pos = Point3D(
                    (line.start.x + line.end.x) / 2,
                    (line.start.y + line.end.y) / 2,
                    (line.start.z + line.end.z) / 2,
                )

                # First part of bond.
                if look == Look.CARTOON:
                    depiction: Depiction = Cartoon(
                        start_color, stroke_color, stroke_width, bond.opacity
                    )

                elif look == Look.GLOSSY:
                    depiction = Glossy(start_color, bond.opacity)

                else:
                    raise ValueError(f"Unknown look: '{look}'")

                scene.add_node(
                    ModelCylinder(
                        Cylinder(line.start, middel_pos, temp_bond_radius, cap_type), depiction
                    )
                )

                # Second part of bond.
                if look == Look.CARTOON:
                    depiction = Cartoon(end_color, stroke_color, stroke_width, bond.opacity)

                elif look == Look.GLOSSY:
                    depiction = Glossy(end_color, bond.opacity)

                else:
                    raise ValueError(f"Unknown look: '{look}'")

                scene.add_node(
                    ModelCylinder(
                        Cylinder(middel_pos, line.end, temp_bond_radius, cap_type), depiction
                    )
                )

            else:
                # Bond as a whole.
                if look == Look.CARTOON:
                    depiction = Cartoon(start_color, stroke_color, stroke_width, bond.opacity)

                elif look == Look.GLOSSY:
                    depiction = Glossy(start_color, bond.opacity)

                else:
                    raise ValueError(f"Unknown look: '{look}'")

                scene.add_node(
                    ModelCylinder(
                        Cylinder(line.start, line.end, temp_bond_radius, cap_type), depiction
                    )
                )


def draw_molecule(
    atoms: ty.List[Atom],
    bonds: ty.List[Bond],
    style: Style,
    look: Look,
    resolution: int,
    window: ty.Optional[ty.Tuple[float, float]] = None,
    view_box: ty.Optional[ty.Tuple[float, float, float, float]] = None,
    rotation_over_x_axis: float = 0.0,
    rotation_over_y_axis: float = 0.0,
    rotation_over_z_axis: float = 0.0,
    scale: float = 1.0,
    focal_length: ty.Optional[float] = None,
    exclude_atoms: ty.Optional[ty.List[str]] = None,
) -> Svg:
    """
    Draw a molecule using the given atoms and bonds.

    :param atoms: The atoms in the molecule.
    :type atoms: ty.List[Atom]
    :param bonds: The bonds in the molecule.
    :type bonds: ty.List[Bond]
    :param style: The style of the depiction.
    :type style: Style
    :param look: The look of the depiction.
    :type look: Look
    :param resolution: The resolution of the depiction.
    :type resolution: int
    :param window: The window of the depiction.
    :type window: ty.Optional[ty.Tuple[float, float]]
    :param view_box: The view box of the depiction.
    :type view_box: ty.Optional[ty.Tuple[float, float, float, float]]
    :param rotation_over_x_axis: The rotation over the x-axis of the depiction.
    :type rotation_over_x_axis: float
    :param rotation_over_y_axis: The rotation over the y-axis of the depiction.
    :type rotation_over_y_axis: float
    :param rotation_over_z_axis: The rotation over the z-axis of the depiction.
    :type rotation_over_z_axis: float
    :param scale: The scale of the depiction.
    :param focal_length: The focal length of the depiction. If None, the focal length
        is calculated based on the dimensions of the scene.
    :type focal_length: ty.Optional[float]
    :param exclude_atoms: The atoms to exclude from the depiction (list of atom symbols).
    :type exclude_atoms: ty.Optional[ty.List[str]]
    :return: The SVG of the depiction.
    :rtype: Svg
    :raises ValueError: If an unknown style or look is given.
    """
    atoms, bonds = filter_atoms_and_bonds(atoms, bonds, exclude_atoms)

    scene = Scene(nodes=[])  # Define empty scene explicitly to prevent memory leak.

    # Default settings for drawing.
    include_spheres = False
    include_cylinders = False
    include_wires = False
    calculate_sphere_sphere_intersections = False
    calculate_sphere_cylinder_intersections = False
    calculate_cylinder_sphere_intersections = False
    calculate_cylinder_cylinder_intersections = False

    # Default settings for cartoon look.
    stroke_color = Color(0, 0, 0)
    stroke_width = 0.05

    # Wire style has a separate implementation that is faster than for geometric shapes.
    if style == Style.WIREFRAME:
        include_wires = True
        draw_bonds_in_wireframe_style(scene, atoms, bonds)

    elif style == Style.SPACEFILLING:
        include_spheres = True
        calculate_sphere_sphere_intersections = True
        draw_atoms_in_spacefilling_style(scene, atoms, look, stroke_color, stroke_width)

    elif style == Style.BALL_AND_STICK or style == Style.TUBE:

        if style == Style.BALL_AND_STICK:
            # We need to draw atoms a sspheres and bonds as cylinders for ball-and-stick
            # style.
            include_spheres = True
            include_cylinders = True
            calculate_cylinder_sphere_intersections = True
            cap_type = CylinderCapType.NO_CAP
            draw_atoms_in_spacefilling_style(
                scene,
                atoms,
                look,
                stroke_color,
                stroke_width,
                radius_scale=1.0 / 3.0,  # Lower the radius of the atoms for ball-and-stick style.
            )

        elif style == Style.TUBE:
            # For tube we only draw bonds as cylinders and therefore the intersections
            # between spheres and cylinders are not needed. However, we do need to
            # calculate the intersections between cylinders.
            include_cylinders = True
            calculate_cylinder_cylinder_intersections = True
            cap_type = CylinderCapType.ROUND

        draw_bonds_in_tube_style(
            scene, atoms, bonds, style, look, cap_type, stroke_color, stroke_width
        )

    else:
        raise ValueError(f"Unknown style: '{style}'")

    # Draw scene.
    svg = scene.draw(
        resolution=resolution,
        window=window,
        view_box=ViewBox(*view_box) if view_box is not None else None,
        rotation_over_x_axis=rotation_over_x_axis,
        rotation_over_y_axis=rotation_over_y_axis,
        rotation_over_z_axis=rotation_over_z_axis,
        include_spheres=include_spheres,
        include_cylinders=include_cylinders,
        include_wires=include_wires,
        calculate_sphere_sphere_intersections=calculate_sphere_sphere_intersections,
        calculate_sphere_cylinder_intersections=calculate_sphere_cylinder_intersections,
        calculate_cylinder_sphere_intersections=calculate_cylinder_sphere_intersections,
        calculate_cylinder_cylinder_intersections=calculate_cylinder_cylinder_intersections,
        filter_nodes_for_intersecting=True,
        scale=scale,
        focal_length=focal_length,
    )

    return svg
