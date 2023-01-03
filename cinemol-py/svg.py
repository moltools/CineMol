from __future__ import annotations
from math import pow
from typing import (Tuple, Callable, Optional, Any, List)
from .fable_modules.fable_library.array import map as map_1
from .fable_modules.fable_library.double import (divide, sqrt)
from .fable_modules.fable_library.list import (map, FSharpList, is_empty, head, tail, filter, empty, append as append_1, singleton as singleton_1, contains)
from .fable_modules.fable_library.reflection import (TypeInfo, union_type)
from .fable_modules.fable_library.seq import (to_list, delay, collect, append, singleton, empty as empty_1)
from .fable_modules.fable_library.string import join
from .fable_modules.fable_library.system_text import (StringBuilder__Append_Z721C83C5, StringBuilder__ctor)
from .fable_modules.fable_library.types import (Array, Union, to_string)
from .fable_modules.fable_library.util import (number_hash, IEnumerable_1, to_enumerable)
from .geometry import (calc_slope, same_side_of_line)
from .helpers import (float_to_str, round as round_1, int_to_str)
from .styles import (Color__Diffuse_5E38073B, atom_color_gradient, get_atom_color, AtomColorStyle)
from .types import (ProjectedAtomInfo, ClipPath, Point2D, Point2D__FindVector_591E284C, SelectForSide, Point2D__Midpoint_591E284C, Point2D__Distance_591E284C, BondInfo, BondType, Depiction, ProjectedMolecule)

def header(_arg1_: float, _arg1__1: float, _arg1__2: float, _arg1__3: float) -> str:
    _arg: Tuple[float, float, float, float] = (_arg1_, _arg1__1, _arg1__2, _arg1__3)
    return ((((((("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg\n\tid=\"Layer_1\"\n\txmlns=\"http://www.w3.org/2000/svg\"\n\tviewBox=\"" + str(_arg[0] * 1.0)) + " ") + str(_arg[1] * 1.0)) + " ") + str(_arg[2] * 1.0)) + " ") + str(_arg[3] * 1.0)) + "\"\n>"


def clipping_to_mask(a: ProjectedAtomInfo, c: ClipPath) -> str:
    def str_1(f: float, a: ProjectedAtomInfo=a, c: ClipPath=c) -> str:
        return float_to_str(f)

    pattern_input: Tuple[Point2D, Point2D] = c.Line
    og_l2: Point2D = pattern_input[1]
    og_l1: Point2D = pattern_input[0]
    pattern_input_1: Tuple[Point2D, Point2D] = (Point2D(og_l1.X - (a.Radius * Point2D__FindVector_591E284C(og_l1, og_l2).X), og_l1.Y - (a.Radius * Point2D__FindVector_591E284C(og_l1, og_l2).Y)), Point2D(og_l2.X - (a.Radius * Point2D__FindVector_591E284C(og_l2, og_l1).X), og_l2.Y - (a.Radius * Point2D__FindVector_591E284C(og_l2, og_l1).Y)))
    l2: Point2D = pattern_input_1[1]
    l1: Point2D = pattern_input_1[0]
    width: float = 2.0 * a.Radius
    perp_slope: float = -1.0 * divide(1.0, calc_slope(l1, l2))
    t: float = divide(width, sqrt(1.0 + pow(perp_slope, 2.0)))
    l1a: Point2D = Point2D(l1.X + t, l1.Y + (perp_slope * t))
    l1b: Point2D = Point2D(l1.X - t, l1.Y - (perp_slope * t))
    l2a: Point2D = Point2D(l2.X + t, l2.Y + (perp_slope * t))
    l2b: Point2D = Point2D(l2.X - t, l2.Y - (perp_slope * t))
    pattern_input_2: Tuple[Optional[Point2D], Optional[Point2D]]
    match_value: SelectForSide = c.SelectForSide
    if match_value.tag == 1:
        pattern_input_2 = ((l1a, l2a)) if same_side_of_line(l1, l2, match_value.fields[0], l1a) else ((l1b, l2b))

    elif match_value.tag == 2:
        pattern_input_2 = (None, None)

    else: 
        pattern_input_2 = ((l1b, l2b)) if same_side_of_line(l1, l2, match_value.fields[0], l1a) else ((l1a, l2a))

    l2new_parallel: Optional[Point2D] = pattern_input_2[1]
    l1new_parallel: Optional[Point2D] = pattern_input_2[0]
    (pattern_matching_result, l1new_parallel_1, l2new_parallel_1) = (None, None, None)
    if l1new_parallel is not None:
        if l2new_parallel is not None:
            pattern_matching_result = 0
            l1new_parallel_1 = l1new_parallel
            l2new_parallel_1 = l2new_parallel

        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        m: Point2D = Point2D__Midpoint_591E284C(og_l1, og_l2)
        t_1: float = divide(0.1 * Point2D__Distance_591E284C(og_l1, og_l2), sqrt(1.0 + pow(perp_slope, 2.0)))
        qa: Point2D = Point2D(m.X + t_1, m.Y + (perp_slope * t_1))
        qb: Point2D = Point2D(m.X - t_1, m.Y - (perp_slope * t_1))
        q: Optional[Point2D]
        match_value_4: SelectForSide = c.SelectForSide
        if match_value_4.tag == 1:
            q = qb if same_side_of_line(og_l1, og_l2, match_value_4.fields[0], qa) else qa

        elif match_value_4.tag == 2:
            q = None

        else: 
            q = qa if same_side_of_line(og_l1, og_l2, match_value_4.fields[0], qa) else qb

        if q is None:
            return ""

        else: 
            q_1: Point2D = q
            return ((((((((((((((((((("<path d=\"\n                M " + str_1(l1.X)) + " ") + str_1(l1.Y)) + "\n                Q ") + str_1(q_1.X)) + " ") + str_1(q_1.Y)) + " ") + str_1(l2.X)) + " ") + str_1(l2.Y)) + "\n                L ") + str_1(l2new_parallel_1.X)) + " ") + str_1(l2new_parallel_1.Y)) + "\n                L ") + str_1(l1new_parallel_1.X)) + " ") + str_1(l1new_parallel_1.Y)) + "\n                Z\"\n            fill=\"black\"/>"


    elif pattern_matching_result == 1:
        return ""



def write_atom_defs(view_box_: float, view_box__1: float, view_box__2: float, view_box__3: float, ball_and_stick: bool, atom: ProjectedAtomInfo) -> str:
    view_box: Tuple[float, float, float, float] = (view_box_, view_box__1, view_box__2, view_box__3)
    pattern_input: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), atom.AtomType), atom_color_gradient[0])
    r1: float = pattern_input[0]
    g1: float = pattern_input[1]
    b1: float = pattern_input[2]
    pattern_input_1: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), atom.AtomType), atom_color_gradient[1])
    r2: float = pattern_input_1[0]
    g2: float = pattern_input_1[1]
    b2: float = pattern_input_1[2]
    pattern_input_2: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), atom.AtomType), atom_color_gradient[2])
    r3: float = pattern_input_2[0]
    g3: float = pattern_input_2[1]
    b3: float = pattern_input_2[2]
    pattern_input_3: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), atom.AtomType), atom_color_gradient[3])
    r4: float = pattern_input_3[0]
    g4: float = pattern_input_3[1]
    b4: float = pattern_input_3[2]
    pattern_input_4: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), atom.AtomType), atom_color_gradient[4])
    r5: float = pattern_input_4[0]
    g5: float = pattern_input_4[1]
    b5: float = pattern_input_4[2]
    if ball_and_stick:
        return ((((((((((((((((((((((((((((((((((((((((((((((((((("\n<radialGradient\n\tid=\"radial-gradient-" + str(atom.Index)) + "\"\n\tcx=\"") + float_to_str(atom.Center.X)) + "\"\n\tcy=\"") + float_to_str(atom.Center.Y)) + "\"\n\tfx=\"") + float_to_str(atom.Center.X)) + "\"\n\tfy=\"") + float_to_str(atom.Center.Y)) + "\"\n\tr=\"") + float_to_str(divide(atom.Radius, 2.0) + 0.4)) + "\"\n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\n\tgradientUnits=\"userSpaceOnUse\"\n>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[0])) + "\" stop-color=\"rgb(") + float_to_str(r1)) + ",") + float_to_str(g1)) + ",") + float_to_str(b1)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[1])) + "\" stop-color=\"rgb(") + float_to_str(r2)) + ",") + float_to_str(g2)) + ",") + float_to_str(b2)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[2])) + "\" stop-color=\"rgb(") + float_to_str(r3)) + ",") + float_to_str(g3)) + ",") + float_to_str(b3)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[3])) + "\" stop-color=\"rgb(") + float_to_str(r4)) + ",") + float_to_str(g4)) + ",") + float_to_str(b4)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[4])) + "\" stop-color=\"rgb(") + float_to_str(r5)) + ",") + float_to_str(g5)) + ",") + float_to_str(b5)) + ")\"/>\n</radialGradient>"

    else: 
        def str_1(f_5: float, view_box_: float=view_box_, view_box__1: float=view_box__1, view_box__2: float=view_box__2, view_box__3: float=view_box__3, ball_and_stick: bool=ball_and_stick, atom: ProjectedAtomInfo=atom) -> str:
            return float_to_str(f_5)

        cx: str = str_1(round_1(3, atom.Center.X))
        cy: str = str_1(round_1(3, atom.Center.Y))
        str_1(round_1(3, atom.Radius))
        def mapping(c: ClipPath, view_box_: float=view_box_, view_box__1: float=view_box__1, view_box__2: float=view_box__2, view_box__3: float=view_box__3, ball_and_stick: bool=ball_and_stick, atom: ProjectedAtomInfo=atom) -> str:
            return clipping_to_mask(atom, c)

        masks: str = join(" ", map(mapping, atom.ClipPaths))
        return (((((((((((((((((((((((((((((((((((((((((((((((((((("\n<radialGradient\n\tid=\"radial-gradient-" + str(atom.Index)) + "\"\n\tcx=\"") + cx) + "\"\n\tcy=\"") + cy) + "\"\n\tfx=\"") + cx) + "\"\n\tfy=\"") + cy) + "\"\n\tr=\"") + float_to_str(atom.Radius + 0.8)) + "\"\n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\n\tgradientUnits=\"userSpaceOnUse\"\n>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[0])) + "\" stop-color=\"rgb(") + float_to_str(r1)) + ",") + float_to_str(g1)) + ",") + float_to_str(b1)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[1])) + "\" stop-color=\"rgb(") + float_to_str(r2)) + ",") + float_to_str(g2)) + ",") + float_to_str(b2)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[2])) + "\" stop-color=\"rgb(") + float_to_str(r3)) + ",") + float_to_str(g3)) + ",") + float_to_str(b3)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[3])) + "\" stop-color=\"rgb(") + float_to_str(r4)) + ",") + float_to_str(g4)) + ",") + float_to_str(b4)) + ")\"/>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[4])) + "\" stop-color=\"rgb(") + float_to_str(r5)) + ",") + float_to_str(g5)) + ",") + float_to_str(b5)) + ")\"/>\n</radialGradient>") + (((((((((((("<mask id=\"mask-" + str(atom.Index)) + "\">\n        <rect id=\"bg\" x=\"") + str_1(view_box[0])) + "\" y=\"") + str_1(view_box[1])) + "\" width=\"") + str_1(view_box[2])) + "\" height=\"") + str_1(view_box[3])) + "\" fill=\"white\"/>\n        ") + masks) + "\n        </mask>")



def write_atoms_defs(vb_: float, vb__1: float, vb__2: float, vb__3: float, atoms: FSharpList[ProjectedAtomInfo], ball_and_stick: bool) -> str:
    vb: Tuple[float, float, float, float] = (vb_, vb__1, vb__2, vb__3)
    def mapping(atom: ProjectedAtomInfo, vb_: float=vb_, vb__1: float=vb__1, vb__2: float=vb__2, vb__3: float=vb__3, atoms: Any=atoms, ball_and_stick: bool=ball_and_stick) -> str:
        return write_atom_defs(vb[0], vb[1], vb[2], vb[3], ball_and_stick, atom)

    return join("", map(mapping, atoms))


def write_atoms_style(atoms: FSharpList[ProjectedAtomInfo]) -> str:
    def mapping(atom: ProjectedAtomInfo, atoms: Any=atoms) -> str:
        return ((("\n.atom-" + str(atom.Index)) + "{fill:url(#radial-gradient-") + str(atom.Index)) + ");}"

    return join("", map(mapping, atoms))


def write_bond_defs(atoms: FSharpList[ProjectedAtomInfo], bond: BondInfo) -> str:
    def find_atom(l_mut: FSharpList[ProjectedAtomInfo], atom_idx_mut: int, atoms: Any=atoms, bond: BondInfo=bond) -> Optional[ProjectedAtomInfo]:
        while True:
            (l, atom_idx) = (l_mut, atom_idx_mut)
            if not is_empty(l):
                x: ProjectedAtomInfo = head(l)
                if x.Index == atom_idx:
                    return x

                else: 
                    l_mut = tail(l)
                    atom_idx_mut = atom_idx
                    continue


            else: 
                return None

            break

    matchValue: Optional[ProjectedAtomInfo] = find_atom(atoms, bond.Start)
    matchValue_1: Optional[ProjectedAtomInfo] = find_atom(atoms, bond.End)
    (pattern_matching_result, e, s) = (None, None, None)
    if matchValue is not None:
        if matchValue_1 is not None:
            pattern_matching_result = 0
            e = matchValue_1
            s = matchValue

        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        pattern_input: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), s.AtomType), atom_color_gradient[2])
        pattern_input_1: Tuple[float, float, float] = Color__Diffuse_5E38073B(get_atom_color(AtomColorStyle(0), e.AtomType), atom_color_gradient[2])
        return ((((((((((((((((((((((((((((((((((((((("\n<linearGradient\n\tid=\"linear-gradient-" + str(bond.Index)) + "-atom-") + str(s.Index)) + "\"\n\tx1=\"") + float_to_str(s.Center.X)) + "\"\n\ty1=\"") + float_to_str(s.Center.Y)) + "\"\n\tx2=\"") + float_to_str(e.Center.X)) + "\"\n\ty2=\"") + float_to_str(e.Center.Y)) + "\"\n\tgradientUnits=\"userSpaceOnUse\"\n>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[2])) + "\" stop-color=\"rgb(") + float_to_str(pattern_input[0])) + ",") + float_to_str(pattern_input[1])) + ",") + float_to_str(pattern_input[2])) + ")\"/>\n</linearGradient>\n        \n<linearGradient\n\tid=\"linear-gradient-") + str(bond.Index)) + "-atom-") + str(e.Index)) + "\"\n\tx1=\"") + float_to_str(s.Center.X)) + "\"\n\ty1=\"") + float_to_str(s.Center.Y)) + "\"\n\tx2=\"") + float_to_str(e.Center.X)) + "\"\n\ty2=\"") + float_to_str(e.Center.Y)) + "\"\n\tgradientUnits=\"userSpaceOnUse\"\n>\n<stop offset=\"") + float_to_str(1.0 - atom_color_gradient[2])) + "\" stop-color=\"rgb(") + float_to_str(pattern_input_1[0])) + ",") + float_to_str(pattern_input_1[1])) + ",") + float_to_str(pattern_input_1[2])) + ")\"/>\n</linearGradient>"

    elif pattern_matching_result == 1:
        return ""



def write_bonds_defs(atoms: FSharpList[ProjectedAtomInfo], bonds: FSharpList[BondInfo]) -> str:
    def mapping(bond: BondInfo, atoms: Any=atoms, bonds: Any=bonds) -> str:
        return write_bond_defs(atoms, bond)

    return join("", map(mapping, bonds))


def write_bonds_style(bonds: FSharpList[BondInfo]) -> str:
    def mapping(bond: BondInfo, bonds: Any=bonds) -> str:
        return ((((((((((((((("\n.bond-" + str(bond.Index)) + "-atom-") + str(bond.Start)) + "{fill:url(#linear-gradient-") + str(bond.Index)) + "-atom-") + str(bond.Start)) + ");}\n        \n.bond-") + str(bond.Index)) + "-atom-") + str(bond.End)) + "{fill:url(#linear-gradient-") + str(bond.Index)) + "-atom-") + str(bond.End)) + ");}"

    return join("", map(mapping, bonds))


def write_atoms_filled(atoms: FSharpList[ProjectedAtomInfo]) -> str:
    def mapping(atom: ProjectedAtomInfo, atoms: Any=atoms) -> str:
        def str_1(f: float, atom: ProjectedAtomInfo=atom) -> str:
            return float_to_str(f)

        pattern_input: Tuple[float, float, float] = (atom.Center.X, atom.Center.Y, atom.Radius)
        return ((((((((("<circle class=\"atom-" + str(atom.Index)) + "\" style=\"\" cx=\"") + str_1(pattern_input[0])) + "\" cy=\"") + str_1(pattern_input[1])) + "\" r=\"") + str_1(pattern_input[2])) + "\" mask=\"url(#mask-") + str(atom.Index)) + ")\"/>"

    return join("", map(mapping, atoms))


def write_atoms_wire(atoms: FSharpList[ProjectedAtomInfo], bonds: FSharpList[BondInfo]) -> str:
    def find_atom(l_mut: FSharpList[ProjectedAtomInfo], atom_idx_mut: int, atoms: Any=atoms, bonds: Any=bonds) -> Optional[ProjectedAtomInfo]:
        while True:
            (l, atom_idx) = (l_mut, atom_idx_mut)
            if not is_empty(l):
                x: ProjectedAtomInfo = head(l)
                if x.Index == atom_idx:
                    return x

                else: 
                    l_mut = tail(l)
                    atom_idx_mut = atom_idx
                    continue


            else: 
                return None

            break

    def find_bonds(l_1: FSharpList[BondInfo], atom_idx_1: int, atoms: Any=atoms, bonds: Any=bonds) -> FSharpList[BondInfo]:
        def predicate(b: BondInfo, l_1: Any=l_1, atom_idx_1: int=atom_idx_1) -> bool:
            return b.Start == atom_idx_1

        return filter(predicate, l_1)

    drawn_atoms: FSharpList[int] = empty()
    def draw_atom(atom: ProjectedAtomInfo, atoms: Any=atoms, bonds: Any=bonds) -> str:
        pattern_input: Tuple[int, int, int] = get_atom_color(AtomColorStyle(0), atom.AtomType).RGB
        return ((((((((((((("<circle\n\tclass=\"atom-" + str(atom.Index)) + "\"\n\tstyle=\"fill:rgb(") + int_to_str(pattern_input[0])) + ",") + int_to_str(pattern_input[1])) + ",") + int_to_str(pattern_input[2])) + ")\"\n        \n\tcx=\"") + float_to_str(atom.Center.X)) + "\"\n\tcy=\"") + float_to_str(atom.Center.Y)) + "\"\n\tr=\"") + float_to_str(0.05)) + "\"\n/>"

    def _arrow68(__unit: None=None, atoms: Any=atoms, bonds: Any=bonds) -> IEnumerable_1[str]:
        def _arrow67(start_atom: ProjectedAtomInfo) -> IEnumerable_1[str]:
            def _arrow66(__unit: None=None) -> IEnumerable_1[str]:
                nonlocal drawn_atoms
                drawn_atoms = append_1(drawn_atoms, singleton_1(start_atom.Index))
                match_value: FSharpList[BondInfo] = find_bonds(bonds, start_atom.Index)
                if is_empty(match_value):
                    return empty_1()

                else: 
                    def _arrow65(atom_bond: BondInfo) -> IEnumerable_1[str]:
                        class ObjectExpr63:
                            @property
                            def Equals(self) -> Callable[[int, int], bool]:
                                def _arrow62(x_1: int, y: int) -> bool:
                                    return x_1 == y

                                return _arrow62

                            @property
                            def GetHashCode(self) -> Callable[[int], int]:
                                return number_hash

                        if not contains(atom_bond.End, drawn_atoms, ObjectExpr63()):
                            match_value_1: Optional[ProjectedAtomInfo] = find_atom(atoms, atom_bond.End)
                            if match_value_1 is None:
                                return empty_1()

                            else: 
                                pattern_input_1: Tuple[ProjectedAtomInfo, ProjectedAtomInfo] = (start_atom, match_value_1)
                                s: ProjectedAtomInfo = pattern_input_1[0]
                                e: ProjectedAtomInfo = pattern_input_1[1]
                                m: Point2D = Point2D__Midpoint_591E284C(s.Center, e.Center)
                                pattern_input_2: Tuple[int, int, int] = get_atom_color(AtomColorStyle(0), s.AtomType).RGB
                                pattern_input_3: Tuple[int, int, int] = get_atom_color(AtomColorStyle(0), e.AtomType).RGB
                                def _arrow64(__unit: None=None) -> IEnumerable_1[str]:
                                    return singleton(((((((((((((("<line\n\tx1=\"" + str(m.X)) + "\"\n\tx2=\"") + str(e.Center.X)) + "\"\n\ty1=\"") + str(m.Y)) + "\"\n\ty2=\"") + str(e.Center.Y)) + "\"\n\tstyle=\"stroke:rgb(") + int_to_str(pattern_input_3[0])) + ",") + int_to_str(pattern_input_3[1])) + ",") + int_to_str(pattern_input_3[2])) + ");stroke-width:0.1\"/>")

                                return append(singleton(((((((((((((("<line\n\tx1=\"" + str(s.Center.X)) + "\"\n\tx2=\"") + str(m.X)) + "\"\n\ty1=\"") + str(s.Center.Y)) + "\"\n\ty2=\"") + str(m.Y)) + "\"\n\tstyle=\"stroke:rgb(") + int_to_str(pattern_input_2[0])) + ",") + int_to_str(pattern_input_2[1])) + ",") + int_to_str(pattern_input_2[2])) + ");stroke-width:0.1\"/>"), delay(_arrow64))


                        else: 
                            return empty_1()


                    return collect(_arrow65, match_value)


            return append(singleton(draw_atom(start_atom)), delay(_arrow66))

        return collect(_arrow67, atoms)

    return join("", to_list(delay(_arrow68)))


def _expr69() -> TypeInfo:
    return union_type("Client.CineMol.Svg.BondEnd", [], BondEnd, lambda: [[], []])


class BondEnd(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["Start", "End"]


BondEnd_reflection = _expr69

def write_ball_and_stick(atoms: FSharpList[ProjectedAtomInfo], bonds: FSharpList[BondInfo]) -> str:
    def round(f: float, atoms: Any=atoms, bonds: Any=bonds) -> float:
        return round_1(3, f)

    def draw_atom(atom: ProjectedAtomInfo, atoms: Any=atoms, bonds: Any=bonds) -> str:
        return ((((((("<circle\n\tclass=\"atom-" + str(atom.Index)) + "\"\n\tstyle=\"\"\n        \n\tcx=\"") + float_to_str(round(atom.Center.X))) + "\"\n\tcy=\"") + float_to_str(round(atom.Center.Y))) + "\"\n\tr=\"") + float_to_str(round(divide(atom.Radius, 2.0)))) + "\"\n/>"

    def draw_single_bond(bond: BondInfo, s: ProjectedAtomInfo, e: ProjectedAtomInfo, atoms: Any=atoms, bonds: Any=bonds) -> str:
        width: float = 0.15 * bond.Scaling
        s_proj: Point2D = Point2D(round(s.Center.X), round(s.Center.Y))
        e_proj: Point2D = Point2D(round(e.Center.X), round(e.Center.Y))
        factor: float = 0.2
        s_proj_1: Point2D = Point2D(round(s_proj.X + ((factor * s.Radius) * Point2D__FindVector_591E284C(s_proj, e_proj).X)), round(s_proj.Y + ((factor * s.Radius) * Point2D__FindVector_591E284C(s_proj, e_proj).Y)))
        e_proj_1: Point2D = Point2D(round(e_proj.X + ((factor * e.Radius) * Point2D__FindVector_591E284C(e_proj, s_proj_1).X)), round(e_proj.Y + ((factor * e.Radius) * Point2D__FindVector_591E284C(e_proj, s_proj_1).Y)))
        midpoint: Point2D
        p: Point2D = Point2D__Midpoint_591E284C(s_proj_1, e_proj_1)
        midpoint = Point2D(round(p.X), round(p.Y))
        def mapping(tupled_arg: Tuple[BondEnd, Point2D, Point2D], bond: BondInfo=bond, s: ProjectedAtomInfo=s, e: ProjectedAtomInfo=e) -> str:
            s_proj_2: Point2D = tupled_arg[1]
            e_proj_2: Point2D = tupled_arg[2]
            perp_slope: float = -1.0 * divide(1.0, calc_slope(s_proj_2, e_proj_2))
            t: float = divide(width, sqrt(1.0 + pow(perp_slope, 2.0)))
            s_top: Point2D = Point2D(round(s_proj_2.X + t), round(s_proj_2.Y + (perp_slope * t)))
            s_bot: Point2D = Point2D(round(s_proj_2.X - t), round(s_proj_2.Y - (perp_slope * t)))
            e_top: Point2D = Point2D(round(e_proj_2.X + t), round(e_proj_2.Y + (perp_slope * t)))
            e_bot: Point2D = Point2D(round(e_proj_2.X - t), round(e_proj_2.Y - (perp_slope * t)))
            s_sweep: int = ((1, 0)) if (s_proj_2.Y > e_proj_2.Y) else ((0, 1))[0] or 0
            if tupled_arg[0].tag == 1:
                return ((((((((((((((((((((((((((((("<path \n\tclass=\"bond-" + str(bond.Index)) + "-atom-") + str(e.Index)) + "\"\n                \n\tstyle=\"\"\n                \n\td= \n\t\t\"M ") + float_to_str(s_top.X)) + " ") + float_to_str(s_top.Y)) + "\n\t\tA ") + float_to_str(width)) + " ") + float_to_str(width)) + " 0 1 ") + int_to_str(s_sweep)) + " ") + float_to_str(s_bot.X)) + " ") + float_to_str(s_bot.Y)) + "\n\t\tL ") + float_to_str(e_bot.X)) + " ") + float_to_str(e_bot.Y)) + "\n\t\tL ") + float_to_str(e_top.X)) + " ") + float_to_str(e_top.Y)) + "\n\t\tL ") + float_to_str(s_top.X)) + " ") + float_to_str(s_top.Y)) + "\"\n/>"

            else: 
                return ((((((((((((((((((((((((((((("<path \n\tclass=\"bond-" + str(bond.Index)) + "-atom-") + str(s.Index)) + "\"\n                \n\tstyle=\"\"\n                \n\td= \n\t\t\"M ") + float_to_str(s_top.X)) + " ") + float_to_str(s_top.Y)) + "\n\t\tA ") + float_to_str(width)) + " ") + float_to_str(width)) + " 0 1 ") + int_to_str(s_sweep)) + " ") + float_to_str(s_bot.X)) + " ") + float_to_str(s_bot.Y)) + "\n\t\tL ") + float_to_str(e_bot.X)) + " ") + float_to_str(e_bot.Y)) + "\n\t\tL ") + float_to_str(e_top.X)) + " ") + float_to_str(e_top.Y)) + "\n\t\tL ") + float_to_str(s_top.X)) + " ") + float_to_str(s_top.Y)) + "\"\n/>"


        return join("", map_1(mapping, [(BondEnd(0), s_proj_1, midpoint), (BondEnd(1), midpoint, e_proj_1)], None))

    def draw_double_bond(bond_1: BondInfo, s_1: ProjectedAtomInfo, e_1: ProjectedAtomInfo, atoms: Any=atoms, bonds: Any=bonds) -> str:
        s_proj_3: Point2D = Point2D(round(s_1.Center.X), round(s_1.Center.Y))
        e_proj_3: Point2D = Point2D(round(e_1.Center.X), round(e_1.Center.Y))
        factor_1: float = 0.2
        s_proj_4: Point2D = Point2D(round(s_proj_3.X + ((factor_1 * s_1.Radius) * Point2D__FindVector_591E284C(s_proj_3, e_proj_3).X)), round(s_proj_3.Y + ((factor_1 * s_1.Radius) * Point2D__FindVector_591E284C(s_proj_3, e_proj_3).Y)))
        e_proj_4: Point2D = Point2D(round(e_proj_3.X + ((factor_1 * e_1.Radius) * Point2D__FindVector_591E284C(e_proj_3, s_proj_4).X)), round(e_proj_3.Y + ((factor_1 * e_1.Radius) * Point2D__FindVector_591E284C(e_proj_3, s_proj_4).Y)))
        def mapping_1(tupled_arg_1: Tuple[BondEnd, Point2D, Point2D], bond_1: BondInfo=bond_1, s_1: ProjectedAtomInfo=s_1, e_1: ProjectedAtomInfo=e_1) -> str:
            s_proj_5: Point2D = tupled_arg_1[1]
            e_proj_5: Point2D = tupled_arg_1[2]
            width_1: float = 0.1 * bond_1.Scaling
            perp_slope_1: float = -1.0 * divide(1.0, calc_slope(s_proj_5, e_proj_5))
            t_1: float = divide(width_1, sqrt(1.0 + pow(perp_slope_1, 2.0)))
            pattern_input_1: Tuple[Point2D, Point2D] = (Point2D(round(s_proj_5.X + t_1), round(s_proj_5.Y + (perp_slope_1 * t_1))), Point2D(round(e_proj_5.X + t_1), round(e_proj_5.Y + (perp_slope_1 * t_1))))
            pattern_input_2: Tuple[Point2D, Point2D] = (Point2D(round(s_proj_5.X - t_1), round(s_proj_5.Y - (perp_slope_1 * t_1))), Point2D(round(e_proj_5.X - t_1), round(e_proj_5.Y - (perp_slope_1 * t_1))))
            def construct_cylinder(p1: Point2D, p2: Point2D, tupled_arg_1: Any=tupled_arg_1) -> Tuple[Point2D, Point2D, Point2D, Point2D, int, int]:
                width_2: float = 0.05 * bond_1.Scaling
                perp_slope_2: float = -1.0 * divide(1.0, calc_slope(p1, p2))
                t_2: float = divide(width_2, sqrt(1.0 + pow(perp_slope_2, 2.0)))
                pattern_input_3: Tuple[int, int] = ((1, 0)) if (p1.Y > p2.Y) else ((0, 1))
                return (Point2D(round(p1.X + t_2), round(p1.Y + (perp_slope_2 * t_2))), Point2D(round(p1.X - t_2), round(p1.Y - (perp_slope_2 * t_2))), Point2D(round(p2.X + t_2), round(p2.Y + (perp_slope_2 * t_2))), Point2D(round(p2.X - t_2), round(p2.Y - (perp_slope_2 * t_2))), pattern_input_3[0], pattern_input_3[1])

            pattern_input_4: Tuple[Point2D, Point2D, Point2D, Point2D, int, int] = construct_cylinder(pattern_input_1[0], pattern_input_1[1])
            s_top1: Point2D = pattern_input_4[0]
            s_bot1: Point2D = pattern_input_4[1]
            e_top1: Point2D = pattern_input_4[2]
            e_bot1: Point2D = pattern_input_4[3]
            pattern_input_5: Tuple[Point2D, Point2D, Point2D, Point2D, int, int] = construct_cylinder(pattern_input_2[0], pattern_input_2[1])
            s_top2: Point2D = pattern_input_5[0]
            s_bot2: Point2D = pattern_input_5[1]
            e_top2: Point2D = pattern_input_5[2]
            e_bot2: Point2D = pattern_input_5[3]
            bond_end_4: int = (e_1.Index if (tupled_arg_1[0].tag == 1) else s_1.Index) or 0
            return join("", to_enumerable([((((((((((((((((((((((((((((("<path \n\tclass=\"bond-" + str(bond_1.Index)) + "-atom-") + str(bond_end_4)) + "\"\n                \n\tstyle=\"\"\n                \n\td= \n\t\t\"M ") + float_to_str(s_top1.X)) + " ") + float_to_str(s_top1.Y)) + "\n\t\tA ") + float_to_str(width_1)) + " ") + float_to_str(width_1)) + " 0 0 ") + int_to_str(pattern_input_4[4])) + " ") + float_to_str(s_bot1.X)) + " ") + float_to_str(s_bot1.Y)) + "\n\t\tL ") + float_to_str(e_bot1.X)) + " ") + float_to_str(e_bot1.Y)) + "\n\t\tL ") + float_to_str(e_top1.X)) + " ") + float_to_str(e_top1.Y)) + "\n\t\tL ") + float_to_str(s_top1.X)) + " ") + float_to_str(s_top1.Y)) + "\"\n/>", ((((((((((((((((((((((((((((("<path \n\tclass=\"bond-" + str(bond_1.Index)) + "-atom-") + str(bond_end_4)) + "\"\n                \n\tstyle=\"\"\n                \n\td= \n\t\t\"M ") + float_to_str(s_top2.X)) + " ") + float_to_str(s_top2.Y)) + "\n\t\tA ") + float_to_str(width_1)) + " ") + float_to_str(width_1)) + " 0 0 ") + int_to_str(pattern_input_5[4])) + " ") + float_to_str(s_bot2.X)) + " ") + float_to_str(s_bot2.Y)) + "\n\t\tL ") + float_to_str(e_bot2.X)) + " ") + float_to_str(e_bot2.Y)) + "\n\t\tL ") + float_to_str(e_top2.X)) + " ") + float_to_str(e_top2.Y)) + "\n\t\tL ") + float_to_str(s_top2.X)) + " ") + float_to_str(s_top2.Y)) + "\"\n/>"]))

        return join("", map_1(mapping_1, [(BondEnd(0), s_proj_4, Point2D__Midpoint_591E284C(s_proj_4, e_proj_4)), (BondEnd(1), Point2D__Midpoint_591E284C(s_proj_4, e_proj_4), e_proj_4)], None))

    def find_atom(l_mut: FSharpList[ProjectedAtomInfo], atom_idx_mut: int, atoms: Any=atoms, bonds: Any=bonds) -> Optional[ProjectedAtomInfo]:
        while True:
            (l, atom_idx) = (l_mut, atom_idx_mut)
            if not is_empty(l):
                x: ProjectedAtomInfo = head(l)
                if x.Index == atom_idx:
                    return x

                else: 
                    l_mut = tail(l)
                    atom_idx_mut = atom_idx
                    continue


            else: 
                return None

            break

    def find_bonds(l_1: FSharpList[BondInfo], atom_idx_1: int, atoms: Any=atoms, bonds: Any=bonds) -> FSharpList[BondInfo]:
        def predicate(b: BondInfo, l_1: Any=l_1, atom_idx_1: int=atom_idx_1) -> bool:
            return b.Start == atom_idx_1

        return filter(predicate, l_1)

    drawn_atoms: FSharpList[int] = empty()
    def _arrow75(__unit: None=None, atoms: Any=atoms, bonds: Any=bonds) -> IEnumerable_1[str]:
        def _arrow74(start_atom: ProjectedAtomInfo) -> IEnumerable_1[str]:
            def _arrow73(__unit: None=None) -> IEnumerable_1[str]:
                nonlocal drawn_atoms
                drawn_atoms = append_1(drawn_atoms, singleton_1(start_atom.Index))
                match_value: FSharpList[BondInfo] = find_bonds(bonds, start_atom.Index)
                if is_empty(match_value):
                    return empty_1()

                else: 
                    def _arrow72(atom_bond: BondInfo) -> IEnumerable_1[str]:
                        class ObjectExpr71:
                            @property
                            def Equals(self) -> Callable[[int, int], bool]:
                                def _arrow70(x_1: int, y: int) -> bool:
                                    return x_1 == y

                                return _arrow70

                            @property
                            def GetHashCode(self) -> Callable[[int], int]:
                                return number_hash

                        if not contains(atom_bond.End, drawn_atoms, ObjectExpr71()):
                            match_value_1: Optional[ProjectedAtomInfo] = find_atom(atoms, atom_bond.End)
                            if match_value_1 is None:
                                return empty_1()

                            else: 
                                end_atom: ProjectedAtomInfo = match_value_1
                                match_value_2: BondType = atom_bond.BondType
                                if match_value_2.tag == 0:
                                    return singleton(draw_single_bond(atom_bond, start_atom, end_atom))

                                elif match_value_2.tag == 1:
                                    return singleton(draw_double_bond(atom_bond, start_atom, end_atom))

                                else: 
                                    return empty_1()



                        else: 
                            return empty_1()


                    return collect(_arrow72, match_value)


            return append(singleton(draw_atom(start_atom)), delay(_arrow73))

        return collect(_arrow74, atoms)

    return join("", to_list(delay(_arrow75)))


def add(s: str, sb: Any) -> Any:
    return StringBuilder__Append_Z721C83C5(sb, s)


def stringify(sb: Any) -> str:
    return to_string(sb)


def write_svg(view_box_: float, view_box__1: float, view_box__2: float, view_box__3: float, depiction: Depiction, mol: ProjectedMolecule) -> str:
    view_box: Tuple[float, float, float, float] = (view_box_, view_box__1, view_box__2, view_box__3)
    if depiction.tag == 1:
        def _arrow79(__unit: None=None, view_box_: float=view_box_, view_box__1: float=view_box__1, view_box__2: float=view_box__2, view_box__3: float=view_box__3, depiction: Depiction=depiction, mol: ProjectedMolecule=mol) -> Any:
            def _arrow78(__unit: None=None) -> Any:
                sb_15: Any
                def _arrow77(__unit: None=None) -> Any:
                    sb_12: Any
                    def _arrow76(__unit: None=None) -> Any:
                        sb_9: Any = StringBuilder__ctor()
                        return add(header(view_box[0], view_box[1], view_box[2], view_box[3]), sb_9)

                    sb_11: Any = add("\n<defs>\n<style>", _arrow76())
                    sb_12 = add(write_atoms_style(mol.Atoms), sb_11)
                    return add(write_bonds_style(mol.Bonds), sb_12)

                sb_14: Any = add("\n</style>", _arrow77())
                sb_15 = add(write_atoms_defs(view_box[0], view_box[1], view_box[2], view_box[3], mol.Atoms, True), sb_14)
                return add(write_bonds_defs(mol.Atoms, mol.Bonds), sb_15)

            sb_17: Any = add("\n</defs>", _arrow78())
            return add(write_ball_and_stick(mol.Atoms, mol.Bonds), sb_17)

        return stringify(add("\n</svg>", _arrow79()))

    elif depiction.tag == 2:
        def _arrow81(__unit: None=None, view_box_: float=view_box_, view_box__1: float=view_box__1, view_box__2: float=view_box__2, view_box__3: float=view_box__3, depiction: Depiction=depiction, mol: ProjectedMolecule=mol) -> Any:
            def _arrow80(__unit: None=None) -> Any:
                sb_20: Any = StringBuilder__ctor()
                return add(header(view_box[0], view_box[1], view_box[2], view_box[3]), sb_20)

            sb_24: Any = add("\n</defs>", add("\n</style>", add("\n<defs>\n<style>", _arrow80())))
            return add(write_atoms_wire(mol.Atoms, mol.Bonds), sb_24)

        return stringify(add("\n</svg>", _arrow81()))

    else: 
        def _arrow85(__unit: None=None, view_box_: float=view_box_, view_box__1: float=view_box__1, view_box__2: float=view_box__2, view_box__3: float=view_box__3, depiction: Depiction=depiction, mol: ProjectedMolecule=mol) -> Any:
            def _arrow84(__unit: None=None) -> Any:
                def _arrow83(__unit: None=None) -> Any:
                    def _arrow82(__unit: None=None) -> Any:
                        sb: Any = StringBuilder__ctor()
                        return add(header(view_box[0], view_box[1], view_box[2], view_box[3]), sb)

                    sb_2: Any = add("\n<defs>\n<style>", _arrow82())
                    return add(write_atoms_style(mol.Atoms), sb_2)

                sb_4: Any = add("\n</style>", _arrow83())
                return add(write_atoms_defs(view_box[0], view_box[1], view_box[2], view_box[3], mol.Atoms, False), sb_4)

            sb_6: Any = add("\n</defs>", _arrow84())
            return add(write_atoms_filled(mol.Atoms), sb_6)

        return stringify(add("\n</svg>", _arrow85()))



__all__ = ["header", "clipping_to_mask", "write_atom_defs", "write_atoms_defs", "write_atoms_style", "write_bond_defs", "write_bonds_defs", "write_bonds_style", "write_atoms_filled", "write_atoms_wire", "BondEnd_reflection", "write_ball_and_stick", "add", "stringify", "write_svg"]

