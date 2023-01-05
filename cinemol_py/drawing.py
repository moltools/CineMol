from __future__ import annotations
from dataclasses import dataclass
from typing import (Any, Callable, Tuple, Optional)
from .fable_modules.fable_library.double import divide
from .fable_modules.fable_library.list import (FSharpList, filter, map, max, sort_by, empty, zip, append, singleton)
from .fable_modules.fable_library.reflection import (TypeInfo, bool_type, record_type)
from .fable_modules.fable_library.types import Record
from .fable_modules.fable_library.util import (equals, compare_primitives)
from .geometry import clip
from .helpers import round
from .styles import AtomType
from .svg import write_svg
from .types import (Depiction, Depiction_reflection, AtomInfo, Axis, AtomInfo__Rotate, Molecule, Point3D, BondInfo, Rotation, Zoom, Point3D__Distance_591E286D, origin, Vector3D, Vector3D__Cross_4C306638, project as project_1, Point2D, ProjectedAtomInfo, ProjectedMolecule)

def _expr56() -> TypeInfo:
    return record_type("Client.CineMol.Drawing.DrawOptions", [], DrawOptions, lambda: [("Depiction", Depiction_reflection()), ("ShowHydrogenAtoms", bool_type)])


@dataclass(eq = False, repr = False)
class DrawOptions(Record):
    Depiction: Depiction
    ShowHydrogenAtoms: bool

DrawOptions_reflection = _expr56

def DrawOptions_get_init(__unit: None=None) -> DrawOptions:
    return DrawOptions(Depiction(0), False)


def filter_atoms(atom_type: AtomType, atoms: FSharpList[AtomInfo]) -> FSharpList[AtomInfo]:
    def _arrow57(atom: AtomInfo, atom_type: AtomType=atom_type, atoms: Any=atoms) -> bool:
        return not equals(atom.AtomType, atom_type)

    return filter(_arrow57, atoms)


def rotate_atoms(axis: Axis, rads: float, atoms: FSharpList[AtomInfo]) -> FSharpList[AtomInfo]:
    def mapping(atom: AtomInfo, axis: Axis=axis, rads: float=rads, atoms: Any=atoms) -> AtomInfo:
        return AtomInfo__Rotate(atom, axis, (divide(rads, 100.0) * 2.0) * 3.141592653589793)

    return map(mapping, atoms)


def change_distance_to_atoms(ratio: float, mol: Molecule) -> Molecule:
    def transform_atom(a: AtomInfo, ratio: float=ratio, mol: Molecule=mol) -> AtomInfo:
        return AtomInfo(a.Index, a.AtomType, Point3D(a.Center.X * ratio, a.Center.Y * ratio, a.Center.Z * ratio), a.Radius * ratio, a.Color)

    def transform_bond(b: BondInfo, ratio: float=ratio, mol: Molecule=mol) -> BondInfo:
        return BondInfo(b.Index, b.Start, b.End, b.BondType, b.Scaling * ratio, b.StartColor, b.EndColor)

    def _arrow58(a_1: AtomInfo, ratio: float=ratio, mol: Molecule=mol) -> AtomInfo:
        return transform_atom(a_1)

    def _arrow59(b_1: BondInfo, ratio: float=ratio, mol: Molecule=mol) -> BondInfo:
        return transform_bond(b_1)

    return Molecule(map(_arrow58, mol.Atoms), map(_arrow59, mol.Bonds))


def draw(view_box: Optional[Tuple[float, float, float, float]], options: DrawOptions, rotation: Rotation, zoom: Zoom, mol: Molecule) -> Tuple[str, Tuple[float, float, float, float]]:
    mol_1: Molecule = Molecule(rotate_atoms(Axis(0), rotation.AxisY, rotate_atoms(Axis(2), rotation.AxisY, rotate_atoms(Axis(1), rotation.AxisX, mol.Atoms))), mol.Bonds)
    offset_view_box: float
    minimum_offset: float = 2.0
    def mapping(a: AtomInfo, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> float:
        return Point3D__Distance_591E286D(a.Center, origin)

    class ObjectExpr60:
        @property
        def Compare(self) -> Callable[[float, float], int]:
            return compare_primitives

    match_value: float = round(0, 2.0 * max(map(mapping, mol_1.Atoms), ObjectExpr60()))
    offset_view_box = minimum_offset if (match_value < minimum_offset) else match_value
    view_box_1: Tuple[float, float, float, float] = view_box if (view_box is not None) else ((-offset_view_box, -offset_view_box, offset_view_box * 2.0, offset_view_box * 2.0))
    focal_length: float = offset_view_box
    pov: Point3D = Point3D(1E-05, 1E-05, focal_length)
    
    
    dist_pov_origin: float = Point3D__Distance_591E286D(pov, origin)
    mol_2: Molecule = Molecule(filter_atoms(AtomType(0), mol_1.Atoms) if (options.ShowHydrogenAtoms == False) else mol_1.Atoms, mol_1.Bonds)
    
    
    def projection(atom: AtomInfo, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> float:
        return -1.0 * Point3D__Distance_591E286D(atom.Center, pov)

    class ObjectExpr61:
        @property
        def Compare(self) -> Callable[[float, float], int]:
            return compare_primitives
    
    # added code snippet to sort atoms; wasn't compiled correctly by fable-py
    sorted_atoms = []
    for a in mol_2.Atoms:
        dist = -1.0 * Point3D__Distance_591E286D(a.Center, pov)
        sorted_atoms.append((dist, a))
    sorted_atoms = sorted(sorted_atoms, key=lambda x: x[0])
    sorted_atoms = [a for _, a in sorted_atoms]
    atoms = empty()
    for a in sorted_atoms:
        atoms = append(atoms, singleton(a))

    # mol_3: Molecule = Molecule(sort_by(projection, mol_2.Atoms, ObjectExpr61()), mol_2.Bonds)
    mol_3: Molecule = Molecule(atoms, mol_2.Bonds)
    
    perspective_mol: Molecule = change_distance_to_atoms(zoom.Ratio, mol_3)
    def mapping_1(atom_1: AtomInfo, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> AtomInfo:
        return AtomInfo(atom_1.Index, atom_1.AtomType, atom_1.Center, divide(dist_pov_origin, Point3D__Distance_591E286D(pov, atom_1.Center)) * atom_1.Radius, atom_1.Color)

    perspective_mol_1: Molecule = Molecule(map(mapping_1, perspective_mol.Atoms), perspective_mol.Bonds)
    camera_forward: Vector3D = Vector3D(-pov.X, -pov.Y, -pov.Z)
    camera_perpendicular: Vector3D = Vector3D(camera_forward.Y, -camera_forward.X, 0.0)
    camera_horizon: Vector3D = Vector3D__Cross_4C306638(camera_forward, camera_perpendicular)
    def project(p: Point3D, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> Point3D:
        return project_1(camera_perpendicular, camera_horizon, camera_forward, pov, focal_length, p)

    def set_perspective_atom(atom_2: AtomInfo, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> AtomInfo:
        return AtomInfo(atom_2.Index, atom_2.AtomType, project(atom_2.Center), atom_2.Radius, atom_2.Color)

    def mapping_2(atom_3: AtomInfo, view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> AtomInfo:
        return set_perspective_atom(atom_3)

    perspective_mol_2: Molecule = Molecule(map(mapping_2, perspective_mol_1.Atoms), perspective_mol_1.Bonds)
    def mapping_3(tupled_arg: Tuple[AtomInfo, AtomInfo], view_box: Any=view_box, options: DrawOptions=options, rotation: Rotation=rotation, zoom: Zoom=zoom, mol: Molecule=mol) -> ProjectedAtomInfo:
        perspective_atom: AtomInfo = tupled_arg[0]
        return ProjectedAtomInfo(perspective_atom.Index, perspective_atom.AtomType, Point2D(perspective_atom.Center.X, perspective_atom.Center.Y), perspective_atom.Radius, clip(pov, perspective_atom, perspective_mol_2, tupled_arg[1], mol_3) if (options.Depiction.tag == 0) else empty(), perspective_atom.Color)

    return (write_svg(view_box_1[0], view_box_1[1], view_box_1[2], view_box_1[3], options.Depiction, ProjectedMolecule(map(mapping_3, zip(perspective_mol_2.Atoms, mol_3.Atoms)), perspective_mol_2.Bonds)), view_box_1)


__all__ = ["DrawOptions_reflection", "DrawOptions_get_init", "filter_atoms", "rotate_atoms", "change_distance_to_atoms", "draw"]

