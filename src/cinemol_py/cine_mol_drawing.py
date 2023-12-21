from __future__ import annotations
from collections.abc import Callable
from math import pow
from typing import Any
from .cine_mol_types import (Chem_Atom, Geometry_Point3D__Rotate, Chem_Molecule, Geometry_Axis, Geometry_Point2D, Geometry_Point2D__CreateVector_2282202F, Geometry_Vector2D, Geometry_Vector2D__Cross_66BEF79A, Geometry_Sphere, Geometry_Point3D, Geometry_Point3D__Dist_2282200E, Geometry_Sphere__PointsOnSphere_Z524259A4, Geometry_Sphere__Intersects_Z460483F4, Geometry_Circle, Svg_Shape, Geometry_Point3D_op_Addition_504401C0, Geometry_Point3D__Mul_5E38073B, Geometry_Point3D_op_Subtraction_504401C0, Geometry_Point3D__Midpoint_2282200E, Color, Chem_BondType, Chem_Bond, Geometry_Axis_Origin, Svg_ViewBox, Chem_AtomType, ModelStyle, Chem_Molecule__GetAtom_Z524259A4, Chem_Molecule__GetBonds, Svg_Header_New, Svg_SVG, Drawing_DrawingOptions)
from .fable_modules.fable_library.double import (sqrt, divide)
from .fable_modules.fable_library.list import (map, min_by, map_indexed, max_by, FSharpList, item, filter, is_empty, of_array, append, reverse, tail, length as length_1, for_all, head, max, sort_by, concat, empty, contains, cons)
from .fable_modules.fable_library.range import range_big_int
from .fable_modules.fable_library.seq import (to_list, delay, map as map_1)
from .fable_modules.fable_library.util import (compare_primitives, IEnumerable_1, equals, number_hash)

def rotate(mol: Chem_Molecule, axis: Geometry_Axis, rad: float) -> Chem_Molecule:
    def _arrow3(atom_1: Chem_Atom, mol: Any=mol, axis: Any=axis, rad: Any=rad) -> Chem_Atom:
        atom: Chem_Atom = atom_1
        return Chem_Atom(atom.Index, atom.Type, atom.Color, atom.Opacity, Geometry_Point3D__Rotate(atom.Position, axis, rad), atom.Radius)

    return Chem_Molecule(map(_arrow3, mol.Atoms), mol.Bonds)


def quick_hull2d(points: FSharpList[Geometry_Point2D]) -> FSharpList[int]:
    def projection(tupled_arg: tuple[int, Geometry_Point2D], points: Any=points) -> float:
        return tupled_arg[1].X

    def mapping(i: int, point: Geometry_Point2D, points: Any=points) -> tuple[int, Geometry_Point2D]:
        return (i, point)

    class ObjectExpr4:
        @property
        def Compare(self) -> Callable[[float, float], int]:
            return compare_primitives

    min_index: int = min_by(projection, map_indexed(mapping, points), ObjectExpr4())[0] or 0
    def projection_1(tupled_arg_1: tuple[int, Geometry_Point2D], points: Any=points) -> float:
        return tupled_arg_1[1].X

    def mapping_1(i_1: int, point_2: Geometry_Point2D, points: Any=points) -> tuple[int, Geometry_Point2D]:
        return (i_1, point_2)

    class ObjectExpr5:
        @property
        def Compare(self) -> Callable[[float, float], int]:
            return compare_primitives

    max_index: int = max_by(projection_1, map_indexed(mapping_1, points), ObjectExpr5())[0] or 0
    def create_range(length: int, points: Any=points) -> FSharpList[int]:
        return to_list(range_big_int(0, 1, length - 1))

    def process_points(tupled_arg_2: tuple[FSharpList[Geometry_Point2D], FSharpList[int], int, int], points: Any=points) -> FSharpList[int]:
        point_list: FSharpList[Geometry_Point2D] = tupled_arg_2[0]
        a: int = tupled_arg_2[2] or 0
        b: int = tupled_arg_2[3] or 0
        W: Geometry_Vector2D = Geometry_Point2D__CreateVector_2282202F(item(b, point_list), item(a, point_list))
        def predicate(tupled_arg_4: tuple[int, float], tupled_arg_2: Any=tupled_arg_2) -> bool:
            i_4: int = tupled_arg_4[0] or 0
            if (i_4 != a) if (tupled_arg_4[1] > 0.0) else False:
                return i_4 != b

            else: 
                return False


        def mapping_3(tupled_arg_3: tuple[int, Geometry_Point2D], tupled_arg_2: Any=tupled_arg_2) -> tuple[int, float]:
            return (tupled_arg_3[0], Geometry_Vector2D__Cross_66BEF79A(Geometry_Point2D__CreateVector_2282202F(tupled_arg_3[1], item(a, point_list)), W))

        def mapping_2(i_2: int, tupled_arg_2: Any=tupled_arg_2) -> tuple[int, Geometry_Point2D]:
            return (i_2, item(i_2, point_list))

        signed_dist: FSharpList[tuple[int, float]] = filter(predicate, map(mapping_3, map(mapping_2, tupled_arg_2[1])))
        if is_empty(signed_dist):
            return of_array([a, b])

        else: 
            def projection_2(tupled_arg_5: tuple[int, float], tupled_arg_2: Any=tupled_arg_2) -> float:
                return tupled_arg_5[1]

            class ObjectExpr6:
                @property
                def Compare(self) -> Callable[[float, float], int]:
                    return compare_primitives

            max_dist_index: int = max_by(projection_2, signed_dist, ObjectExpr6())[0] or 0
            def _arrow7(tuple_3: tuple[int, float], tupled_arg_2: Any=tupled_arg_2) -> int:
                return tuple_3[0]

            new_index_list: FSharpList[int] = map(_arrow7, signed_dist)
            left: FSharpList[int] = process_points((point_list, new_index_list, a, max_dist_index))
            right: FSharpList[int] = process_points((point_list, new_index_list, max_dist_index, b))
            return append(reverse(tail(reverse(left))), right)


    left_1: FSharpList[int] = process_points((points, create_range(length_1(points)), min_index, max_index))
    right_1: FSharpList[int] = process_points((points, create_range(length_1(points)), max_index, min_index))
    return append(reverse(tail(reverse(left_1))), reverse(tail(reverse(right_1))))


def make_polygon(curr_atom: Chem_Atom, intersects_with: FSharpList[Chem_Atom], resolution: int) -> FSharpList[Geometry_Point2D]:
    sphere: Geometry_Sphere = Geometry_Sphere(curr_atom.Position, curr_atom.Radius)
    def mapping(atom: Chem_Atom, curr_atom: Any=curr_atom, intersects_with: Any=intersects_with, resolution: Any=resolution) -> Geometry_Sphere:
        return Geometry_Sphere(atom.Position, atom.Radius)

    other_spheres: FSharpList[Geometry_Sphere] = map(mapping, intersects_with)
    def mapping_1(point_1: Geometry_Point3D, curr_atom: Any=curr_atom, intersects_with: Any=intersects_with, resolution: Any=resolution) -> Geometry_Point2D:
        return Geometry_Point2D(point_1.X, point_1.Y)

    def not_inside_other_sphere(point: Geometry_Point3D, curr_atom: Any=curr_atom, intersects_with: Any=intersects_with, resolution: Any=resolution) -> bool:
        def predicate(other_sphere: Geometry_Sphere, point: Any=point) -> bool:
            return Geometry_Point3D__Dist_2282200E(other_sphere.Center, point) > other_sphere.Radius

        return for_all(predicate, other_spheres)

    points: FSharpList[Geometry_Point2D] = map(mapping_1, filter(not_inside_other_sphere, Geometry_Sphere__PointsOnSphere_Z524259A4(sphere, resolution)))
    def mapping_2(i: int, curr_atom: Any=curr_atom, intersects_with: Any=intersects_with, resolution: Any=resolution) -> Geometry_Point2D:
        return item(i, points)

    return map(mapping_2, quick_hull2d(points))


def atom_to_svg(prev_atoms: FSharpList[Chem_Atom], curr_atom: Chem_Atom, resolution: int, opacity: float) -> Svg_Shape:
    curr_atom_geom: Geometry_Sphere = Geometry_Sphere(curr_atom.Position, curr_atom.Radius)
    def predicate(atom: Chem_Atom, prev_atoms: Any=prev_atoms, curr_atom: Any=curr_atom, resolution: Any=resolution, opacity: Any=opacity) -> bool:
        return Geometry_Sphere__Intersects_Z460483F4(Geometry_Sphere(atom.Position, atom.Radius), curr_atom_geom)

    intersects_with: FSharpList[Chem_Atom] = filter(predicate, prev_atoms)
    if is_empty(intersects_with):
        return Svg_Shape(0, curr_atom.Index, curr_atom.Color, Geometry_Circle(Geometry_Point2D(curr_atom.Position.X, curr_atom.Position.Y), curr_atom.Radius), opacity)

    else: 
        return Svg_Shape(1, curr_atom.Index, curr_atom.Color, make_polygon(curr_atom, intersects_with, resolution), opacity)



def calculate_perpendicular_points(p1: Geometry_Point3D, p2: Geometry_Point3D, radius: float, num_points: int) -> FSharpList[Geometry_Point3D]:
    dx: float = p2.X - p1.X
    dy: float = p2.Y - p1.Y
    len_1: float = sqrt(pow(dx, 2.0) + pow(dy, 2.0))
    ux: float = divide(dx, len_1)
    uy: float = divide(dy, len_1)
    step: float = divide(2.0 * radius, num_points - 1)
    def _arrow9(__unit: None=None, p1: Any=p1, p2: Any=p2, radius: Any=radius, num_points: Any=num_points) -> IEnumerable_1[Geometry_Point3D]:
        def _arrow8(i: int) -> Geometry_Point3D:
            offset: float = (i * step) - radius
            return Geometry_Point3D(p1.X + (offset * uy), p1.Y - (offset * ux), p1.Z)

        return map_1(_arrow8, range_big_int(0, 1, num_points - 1))

    return to_list(delay(_arrow9))


def draw_single_bond(width: float, index1: int, index2: int, color1: Color, color2: Color, p1: Geometry_Point3D, p2: Geometry_Point3D, opacity1: float, opacity2: float) -> FSharpList[Svg_Shape]:
    p1_1: Geometry_Point3D = Geometry_Point3D_op_Addition_504401C0(p1, Geometry_Point3D__Mul_5E38073B(Geometry_Point3D_op_Subtraction_504401C0(p2, p1), 0.3))
    p2_1: Geometry_Point3D = Geometry_Point3D_op_Addition_504401C0(p2, Geometry_Point3D__Mul_5E38073B(Geometry_Point3D_op_Subtraction_504401C0(p1_1, p2), 0.3))
    points1: FSharpList[Geometry_Point3D] = calculate_perpendicular_points(p1_1, p2_1, width, 2)
    mid_points: FSharpList[Geometry_Point3D] = calculate_perpendicular_points(Geometry_Point3D__Midpoint_2282200E(p1_1, p2_1), p2_1, width, 2)
    points2: FSharpList[Geometry_Point3D] = calculate_perpendicular_points(p2_1, p1_1, width, 2)
    pattern_input: tuple[Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D]
    match_value: FSharpList[Geometry_Point3D] = append(points1, reverse(mid_points))
    (pattern_matching_result, b1, b2, b3, b4) = (None, None, None, None, None)
    if not is_empty(match_value):
        if not is_empty(tail(match_value)):
            if not is_empty(tail(tail(match_value))):
                if not is_empty(tail(tail(tail(match_value)))):
                    if is_empty(tail(tail(tail(tail(match_value))))):
                        pattern_matching_result = 0
                        b1 = head(match_value)
                        b2 = head(tail(match_value))
                        b3 = head(tail(tail(match_value)))
                        b4 = head(tail(tail(tail(match_value))))

                    else: 
                        pattern_matching_result = 1


                else: 
                    pattern_matching_result = 1


            else: 
                pattern_matching_result = 1


        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        pattern_input = (Geometry_Point2D(b1.X, b1.Y), Geometry_Point2D(b2.X, b2.Y), Geometry_Point2D(b3.X, b3.Y), Geometry_Point2D(b4.X, b4.Y))

    elif pattern_matching_result == 1:
        pattern_input = (Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0))

    pattern_input_1: tuple[Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D]
    match_value_1: FSharpList[Geometry_Point3D] = append(mid_points, points2)
    (pattern_matching_result_1, b5, b6, b7, b8) = (None, None, None, None, None)
    if not is_empty(match_value_1):
        if not is_empty(tail(match_value_1)):
            if not is_empty(tail(tail(match_value_1))):
                if not is_empty(tail(tail(tail(match_value_1)))):
                    if is_empty(tail(tail(tail(tail(match_value_1))))):
                        pattern_matching_result_1 = 0
                        b5 = head(match_value_1)
                        b6 = head(tail(match_value_1))
                        b7 = head(tail(tail(match_value_1)))
                        b8 = head(tail(tail(tail(match_value_1))))

                    else: 
                        pattern_matching_result_1 = 1


                else: 
                    pattern_matching_result_1 = 1


            else: 
                pattern_matching_result_1 = 1


        else: 
            pattern_matching_result_1 = 1


    else: 
        pattern_matching_result_1 = 1

    if pattern_matching_result_1 == 0:
        pattern_input_1 = (Geometry_Point2D(b5.X, b5.Y), Geometry_Point2D(b6.X, b6.Y), Geometry_Point2D(b7.X, b7.Y), Geometry_Point2D(b8.X, b8.Y))

    elif pattern_matching_result_1 == 1:
        pattern_input_1 = (Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0), Geometry_Point2D(0.0, 0.0))

    def _arrow11(__unit: None=None, width: Any=width, index1: Any=index1, index2: Any=index2, color1: Any=color1, color2: Any=color2, p1: Any=p1, p2: Any=p2, opacity1: Any=opacity1, opacity2: Any=opacity2) -> Svg_Shape:
        tupled_arg: tuple[int, Color, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, float] = (index1, color1, pattern_input[0], pattern_input[1], pattern_input[2], pattern_input[3], opacity1)
        return Svg_Shape(2, tupled_arg[0], tupled_arg[1], tupled_arg[2], tupled_arg[3], tupled_arg[4], tupled_arg[5], tupled_arg[6])

    def _arrow12(__unit: None=None, width: Any=width, index1: Any=index1, index2: Any=index2, color1: Any=color1, color2: Any=color2, p1: Any=p1, p2: Any=p2, opacity1: Any=opacity1, opacity2: Any=opacity2) -> Svg_Shape:
        tupled_arg_1: tuple[int, Color, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, Geometry_Point2D, float] = (index2, color2, pattern_input_1[0], pattern_input_1[1], pattern_input_1[2], pattern_input_1[3], opacity2)
        return Svg_Shape(2, tupled_arg_1[0], tupled_arg_1[1], tupled_arg_1[2], tupled_arg_1[3], tupled_arg_1[4], tupled_arg_1[5], tupled_arg_1[6])

    return of_array([_arrow11(), _arrow12()])


def bond_to_svg(begin_atom: Chem_Atom, end_atom: Chem_Atom, bond: Chem_Bond) -> FSharpList[Svg_Shape]:
    pattern_input: tuple[Chem_Atom, Chem_Atom] = ((end_atom, begin_atom)) if (begin_atom.Position.Z > end_atom.Position.Z) else ((begin_atom, end_atom))
    a2: Chem_Atom = pattern_input[1]
    a1: Chem_Atom = pattern_input[0]
    p1: Geometry_Point3D = a1.Position
    p2: Geometry_Point3D = a2.Position
    pattern_input_1: tuple[int, int, float, float]
    matchValue: int = bond.BeginAtomIndex or 0
    matchValue_1: int = bond.EndAtomIndex or 0
    pattern_input_1 = ((bond.BeginIndex, bond.EndIndex, a1.Opacity, a2.Opacity)) if ((matchValue_1 == a2.Index) if (matchValue == a1.Index) else False) else (((bond.EndIndex, bond.BeginIndex, a2.Opacity, a1.Opacity)) if ((matchValue_1 == a1.Index) if (matchValue == a2.Index) else False) else ((0, 0, 1.0, 1.0)))
    opacity2: float = pattern_input_1[3]
    opacity1: float = pattern_input_1[2]
    index2: int = pattern_input_1[1] or 0
    index1: int = pattern_input_1[0] or 0
    match_value_1: Chem_BondType = bond.Type
    if match_value_1.tag == 3:
        return draw_single_bond(0.2, index1, index2, a1.Color, a2.Color, p1, p2, opacity1, opacity2)

    elif match_value_1.tag == 1:
        radius: float = divide(bond.Radius, 3.0)
        points1: FSharpList[Geometry_Point3D] = calculate_perpendicular_points(p1, p2, radius, 2)
        points2: FSharpList[Geometry_Point3D] = reverse(calculate_perpendicular_points(p2, p1, radius, 2))
        return append(draw_single_bond(0.1, index1, index2, a1.Color, a2.Color, item(0, points1), item(0, points2), opacity1, opacity2), draw_single_bond(0.1, index1, index2, a1.Color, a2.Color, item(1, points1), item(1, points2), opacity1, opacity2))

    elif match_value_1.tag == 2:
        radius_1: float = divide(bond.Radius, 3.0)
        points1_1: FSharpList[Geometry_Point3D] = calculate_perpendicular_points(p1, p2, radius_1, 3)
        points2_1: FSharpList[Geometry_Point3D] = reverse(calculate_perpendicular_points(p2, p1, radius_1, 3))
        return append(draw_single_bond(0.05, index1, index2, a1.Color, a2.Color, item(0, points1_1), item(0, points2_1), opacity1, opacity2), append(draw_single_bond(0.05, index1, index2, a1.Color, a2.Color, item(1, points1_1), item(1, points2_1), opacity1, opacity2), draw_single_bond(0.05, index1, index2, a1.Color, a2.Color, item(2, points1_1), item(2, points2_1), opacity1, opacity2)))

    else: 
        return draw_single_bond(0.2, index1, index2, a1.Color, a2.Color, p1, p2, opacity1, opacity2)



def bond_to_wire(begin_atom: Chem_Atom, end_atom: Chem_Atom) -> FSharpList[Svg_Shape]:
    pattern_input: tuple[Chem_Atom, Chem_Atom] = ((end_atom, begin_atom)) if (begin_atom.Position.Z > end_atom.Position.Z) else ((begin_atom, end_atom))
    a2: Chem_Atom = pattern_input[1]
    a1: Chem_Atom = pattern_input[0]
    p1: Geometry_Point3D = a1.Position
    p2: Geometry_Point3D = a2.Position
    mid_point: Geometry_Point3D = Geometry_Point3D__Midpoint_2282200E(p1, p2)
    mid_point_1: Geometry_Point2D = Geometry_Point2D(mid_point.X, mid_point.Y)
    def _arrow14(__unit: None=None, begin_atom: Any=begin_atom, end_atom: Any=end_atom) -> Svg_Shape:
        tupled_arg: tuple[int, Color, Geometry_Point2D, Geometry_Point2D, float] = (a1.Index, a1.Color, Geometry_Point2D(p1.X, p1.Y), mid_point_1, a1.Opacity)
        return Svg_Shape(3, tupled_arg[0], tupled_arg[1], tupled_arg[2], tupled_arg[3], tupled_arg[4])

    def _arrow15(__unit: None=None, begin_atom: Any=begin_atom, end_atom: Any=end_atom) -> Svg_Shape:
        tupled_arg_1: tuple[int, Color, Geometry_Point2D, Geometry_Point2D, float] = (a2.Index, a2.Color, mid_point_1, Geometry_Point2D(p2.X, p2.Y), a2.Opacity)
        return Svg_Shape(3, tupled_arg_1[0], tupled_arg_1[1], tupled_arg_1[2], tupled_arg_1[3], tupled_arg_1[4])

    return of_array([_arrow14(), _arrow15()])


def draw(mol: Chem_Molecule, options: Drawing_DrawingOptions) -> tuple[Svg_SVG, Drawing_DrawingOptions]:
    origin: Geometry_Point3D = Geometry_Axis_Origin()
    view_box_2: Svg_ViewBox
    def mapping(atom: Chem_Atom, mol: Any=mol, options: Any=options) -> float:
        return Geometry_Point3D__Dist_2282200E(atom.Position, origin)

    class ObjectExpr16:
        @property
        def Compare(self) -> Callable[[float, float], int]:
            return compare_primitives

    offset: float = 5.0 + max(map(mapping, mol.Atoms), ObjectExpr16())
    match_value: Svg_ViewBox | None = options.ViewBox
    view_box_2 = Svg_ViewBox(-offset, -offset, offset * 2.0, offset * 2.0) if (match_value is None) else match_value
    objs: FSharpList[Svg_Shape]
    def predicate(atom_1: Chem_Atom, mol: Any=mol, options: Any=options) -> bool:
        return not equals(atom_1.Type, Chem_AtomType(0))

    atoms: FSharpList[Chem_Atom] = mol.Atoms if options.DisplayHydrogenAtoms else filter(predicate, mol.Atoms)
    match_value_2: ModelStyle = options.ModelStyle
    if match_value_2.tag == 1:
        def mapping_1(atom_4: Chem_Atom, mol: Any=mol, options: Any=options) -> Chem_Atom:
            return Chem_Atom(atom_4.Index, atom_4.Type, atom_4.Color, atom_4.Opacity, atom_4.Position, atom_4.Radius * 0.3)

        def projection_1(atom_3: Chem_Atom, mol: Any=mol, options: Any=options) -> float:
            return atom_3.Position.Z

        class ObjectExpr17:
            @property
            def Compare(self) -> Callable[[float, float], int]:
                return compare_primitives

        atoms_4: FSharpList[Chem_Atom] = map(mapping_1, reverse(sort_by(projection_1, atoms, ObjectExpr17())))
        def process_atoms_1(prev_atoms_2_mut: FSharpList[Chem_Atom], shapes_1_mut: FSharpList[Svg_Shape], mol: Any=mol, options: Any=options) -> FSharpList[Svg_Shape]:
            while True:
                (prev_atoms_2, shapes_1) = (prev_atoms_2_mut, shapes_1_mut)
                if not is_empty(prev_atoms_2):
                    prev_atoms_3: FSharpList[Chem_Atom] = tail(prev_atoms_2)
                    curr_atom_1: Chem_Atom = head(prev_atoms_2)
                    shape_1: Svg_Shape = atom_to_svg(prev_atoms_3, curr_atom_1, options.Resolution, curr_atom_1.Opacity)
                    def mapping_2(atom_5: Chem_Atom, prev_atoms_2: Any=prev_atoms_2, shapes_1: Any=shapes_1) -> int:
                        return atom_5.Index

                    prev_atom_indices: FSharpList[int] = map(mapping_2, prev_atoms_3)
                    prev_atoms_2_mut = prev_atoms_3
                    def mapping_3(bond_1: Chem_Bond, prev_atoms_2: Any=prev_atoms_2, shapes_1: Any=shapes_1) -> FSharpList[Svg_Shape]:
                        begin_atom: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(mol, bond_1.BeginAtomIndex)
                        end_atom: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(mol, bond_1.EndAtomIndex)
                        (pattern_matching_result, b, e) = (None, None, None)
                        if begin_atom is not None:
                            if end_atom is not None:
                                pattern_matching_result = 0
                                b = begin_atom
                                e = end_atom

                            else: 
                                pattern_matching_result = 1


                        else: 
                            pattern_matching_result = 1

                        if pattern_matching_result == 0:
                            return bond_to_svg(Chem_Atom(b.Index, b.Type, b.Color, b.Opacity, b.Position, b.Radius * 0.3), Chem_Atom(e.Index, e.Type, e.Color, e.Opacity, e.Position, e.Radius * 0.3), bond_1)

                        elif pattern_matching_result == 1:
                            return empty()


                    def predicate_1(bond: Chem_Bond, prev_atoms_2: Any=prev_atoms_2, shapes_1: Any=shapes_1) -> bool:
                        class ObjectExpr19:
                            @property
                            def Equals(self) -> Callable[[int, int], bool]:
                                def _arrow18(x_3: int, y_4: int) -> bool:
                                    return x_3 == y_4

                                return _arrow18

                            @property
                            def GetHashCode(self) -> Callable[[int], int]:
                                return number_hash

                        if not contains(bond.BeginAtomIndex, prev_atom_indices, ObjectExpr19()):
                            class ObjectExpr21:
                                @property
                                def Equals(self) -> Callable[[int, int], bool]:
                                    def _arrow20(x_4: int, y_5: int) -> bool:
                                        return x_4 == y_5

                                    return _arrow20

                                @property
                                def GetHashCode(self) -> Callable[[int], int]:
                                    return number_hash

                            return not contains(bond.EndAtomIndex, prev_atom_indices, ObjectExpr21())

                        else: 
                            return False


                    shapes_1_mut = append(cons(shape_1, concat(map(mapping_3, filter(predicate_1, Chem_Molecule__GetBonds(mol, options.DisplayHydrogenAtoms, curr_atom_1.Index))))), shapes_1)
                    continue

                else: 
                    return shapes_1

                break

        process_atoms_1: Callable[[FSharpList[Chem_Atom], FSharpList[Svg_Shape]], FSharpList[Svg_Shape]] = process_atoms_1
        objs = empty() if is_empty(atoms_4) else process_atoms_1(atoms_4, empty())

    elif match_value_2.tag == 2:
        def projection_2(atom_6: Chem_Atom, mol: Any=mol, options: Any=options) -> float:
            return atom_6.Position.Z

        class ObjectExpr22:
            @property
            def Compare(self) -> Callable[[float, float], int]:
                return compare_primitives

        atoms_6: FSharpList[Chem_Atom] = reverse(sort_by(projection_2, atoms, ObjectExpr22()))
        def process_atoms_2(prev_atoms_4_mut: FSharpList[Chem_Atom], shapes_2_mut: FSharpList[Svg_Shape], mol: Any=mol, options: Any=options) -> FSharpList[Svg_Shape]:
            while True:
                (prev_atoms_4, shapes_2) = (prev_atoms_4_mut, shapes_2_mut)
                if not is_empty(prev_atoms_4):
                    prev_atoms_5: FSharpList[Chem_Atom] = tail(prev_atoms_4)
                    def mapping_4(atom_7: Chem_Atom, prev_atoms_4: Any=prev_atoms_4, shapes_2: Any=shapes_2) -> int:
                        return atom_7.Index

                    prev_atom_indices_1: FSharpList[int] = map(mapping_4, prev_atoms_5)
                    prev_atoms_4_mut = prev_atoms_5
                    def mapping_5(bond_3: Chem_Bond, prev_atoms_4: Any=prev_atoms_4, shapes_2: Any=shapes_2) -> FSharpList[Svg_Shape]:
                        begin_atom_1: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(mol, bond_3.BeginAtomIndex)
                        end_atom_1: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(mol, bond_3.EndAtomIndex)
                        (pattern_matching_result_1, b_2, e_2) = (None, None, None)
                        if begin_atom_1 is not None:
                            if end_atom_1 is not None:
                                pattern_matching_result_1 = 0
                                b_2 = begin_atom_1
                                e_2 = end_atom_1

                            else: 
                                pattern_matching_result_1 = 1


                        else: 
                            pattern_matching_result_1 = 1

                        if pattern_matching_result_1 == 0:
                            return bond_to_wire(b_2, e_2)

                        elif pattern_matching_result_1 == 1:
                            return empty()


                    def predicate_2(bond_2: Chem_Bond, prev_atoms_4: Any=prev_atoms_4, shapes_2: Any=shapes_2) -> bool:
                        class ObjectExpr27:
                            @property
                            def Equals(self) -> Callable[[int, int], bool]:
                                def _arrow26(x_6: int, y_7: int) -> bool:
                                    return x_6 == y_7

                                return _arrow26

                            @property
                            def GetHashCode(self) -> Callable[[int], int]:
                                return number_hash

                        if not contains(bond_2.BeginAtomIndex, prev_atom_indices_1, ObjectExpr27()):
                            class ObjectExpr29:
                                @property
                                def Equals(self) -> Callable[[int, int], bool]:
                                    def _arrow28(x_7: int, y_8: int) -> bool:
                                        return x_7 == y_8

                                    return _arrow28

                                @property
                                def GetHashCode(self) -> Callable[[int], int]:
                                    return number_hash

                            return not contains(bond_2.EndAtomIndex, prev_atom_indices_1, ObjectExpr29())

                        else: 
                            return False


                    shapes_2_mut = append(concat(map(mapping_5, filter(predicate_2, Chem_Molecule__GetBonds(mol, options.DisplayHydrogenAtoms, head(prev_atoms_4).Index)))), shapes_2)
                    continue

                else: 
                    return shapes_2

                break

        process_atoms_2: Callable[[FSharpList[Chem_Atom], FSharpList[Svg_Shape]], FSharpList[Svg_Shape]] = process_atoms_2
        objs = empty() if is_empty(atoms_6) else process_atoms_2(atoms_6, empty())

    else: 
        def projection(atom_2: Chem_Atom, mol: Any=mol, options: Any=options) -> float:
            return atom_2.Position.Z

        class ObjectExpr40:
            @property
            def Compare(self) -> Callable[[float, float], int]:
                return compare_primitives

        atoms_1: FSharpList[Chem_Atom] = reverse(sort_by(projection, atoms, ObjectExpr40()))
        def process_atoms(prev_atoms_mut: FSharpList[Chem_Atom], shapes_mut: FSharpList[Svg_Shape], mol: Any=mol, options: Any=options) -> FSharpList[Svg_Shape]:
            processed = 0
            while True:
                (prev_atoms, shapes) = (prev_atoms_mut, shapes_mut)
                if not is_empty(prev_atoms):
                    prev_atoms_1: FSharpList[Chem_Atom] = tail(prev_atoms)
                    curr_atom: Chem_Atom = head(prev_atoms)
                    prev_atoms_mut = prev_atoms_1
                    shapes_mut = cons(atom_to_svg(prev_atoms_1, curr_atom, options.Resolution, curr_atom.Opacity), shapes)
                    processed += 1
                    print(f"{processed}".zfill(5), end="\r")
                    continue

                else: 
                    return shapes

                break


        process_atoms: Callable[[FSharpList[Chem_Atom], FSharpList[Svg_Shape]], FSharpList[Svg_Shape]] = process_atoms
        objs = empty() if is_empty(atoms_1) else process_atoms(atoms_1, empty())

    return (Svg_SVG(Svg_Header_New(), "model", view_box_2, objs, options.ArtStyle), options)


__all__ = ["rotate", "quick_hull2d", "make_polygon", "atom_to_svg", "calculate_perpendicular_points", "draw_single_bond", "bond_to_svg", "bond_to_wire", "draw"]

