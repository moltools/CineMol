from __future__ import annotations
from math import (pow, inf)
from typing import (Tuple, Optional, Callable)
from .fable_modules.fable_library.array import map as map_1
from .fable_modules.fable_library.double import (sqrt, divide)
from .fable_modules.fable_library.list import (FSharpList, map, contains, zip, empty as empty_1)
from .fable_modules.fable_library.seq import (to_list, delay, collect, singleton, empty, to_array)
from .fable_modules.fable_library.types import Array
from .fable_modules.fable_library.util import (IEnumerable_1, number_hash)
from .types import (Point2D, Point3D, AtomInfo, Molecule, Point3D__Distance_591E286D, AtomInfo__Intersects_Z6467FD91, AtomInfo__Intersection_Z6467FD91, SphereSphereIntersection, Point2D__Distance_591E284C, SelectForSide, ClipPath)

def intersection_between_circles(c_p1: Point2D, r_p1: float, c_p2: Point2D, r_p2: float) -> Optional[Tuple[Point2D, Point2D]]:
    d: float = sqrt(pow(c_p2.X - c_p1.X, 2.0) + pow(c_p2.Y - c_p1.Y, 2.0))
    if d > (r_p1 + r_p2):
        return None

    elif (r_p1 == r_p2) if (d == 0.0) else False:
        return None

    elif True if (d < abs(r_p1 - r_p1)) else (d < 1E-05):
        return None

    else: 
        a: float = divide((pow(r_p1, 2.0) - pow(r_p2, 2.0)) + pow(d, 2.0), 2.0 * d)
        h1: float = pow(r_p1, 2.0) - pow(a, 2.0)
        if h1 <= 0.0:
            return None

        else: 
            h2: float = sqrt(h1)
            x2: float = c_p1.X + divide(a * (c_p2.X - c_p1.X), d)
            y2: float = c_p1.Y + divide(a * (c_p2.Y - c_p1.Y), d)
            return (Point2D(x2 + divide(h2 * (c_p2.Y - c_p1.Y), d), y2 - divide(h2 * (c_p2.X - c_p1.X), d)), Point2D(x2 - divide(h2 * (c_p2.Y - c_p1.Y), d), y2 + divide(h2 * (c_p2.X - c_p1.X), d)))




def calc_slope(p1: Point2D, p2: Point2D) -> float:
    return divide(p2.Y - p1.Y, p2.X - p1.X)


def same_side_of_line(line_: Point2D, line__1: Point2D, p1: Point2D, p2: Point2D) -> bool:
    line: Tuple[Point2D, Point2D] = (line_, line__1)
    l2: Point2D = line[1]
    l1: Point2D = line[0]
    def d(p: Point2D, line_: Point2D=line_, line__1: Point2D=line__1, p1: Point2D=p1, p2: Point2D=p2) -> float:
        return ((p.X - l1.X) * (l2.Y - l1.Y)) - ((p.Y - l1.Y) * (l2.X - l1.X))

    if (d(p1) > 0.0) == (d(p2) > 0.0):
        return True

    else: 
        return False



def clip(pov: Point3D, pers_atom: AtomInfo, pers_mol: Molecule, atom: AtomInfo, mol: Molecule) -> FSharpList[ClipPath]:
    dist_pov_atom: float = Point3D__Distance_591E286D(pov, atom.Center)
    def _arrow29(__unit: None=None, pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> IEnumerable_1[AtomInfo]:
        def _arrow28(a: AtomInfo) -> IEnumerable_1[AtomInfo]:
            return singleton(a) if (dist_pov_atom < Point3D__Distance_591E286D(pov, a.Center)) else empty()

        return collect(_arrow28, mol.Atoms)

    atoms_for_clipping: FSharpList[AtomInfo] = to_list(delay(_arrow29))
    def mapping(a_1: AtomInfo, pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> int:
        return a_1.Index

    inds: FSharpList[int] = map(mapping, atoms_for_clipping)
    def _arrow47(__unit: None=None, pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> IEnumerable_1[AtomInfo]:
        def _arrow46(a_2: AtomInfo) -> IEnumerable_1[AtomInfo]:
            class ObjectExpr45:
                @property
                def Equals(self) -> Callable[[int, int], bool]:
                    def _arrow44(x: int, y: int) -> bool:
                        return x == y

                    return _arrow44

                @property
                def GetHashCode(self) -> Callable[[int], int]:
                    return number_hash

            return singleton(a_2) if contains(a_2.Index, inds, ObjectExpr45()) else empty()

        return collect(_arrow46, pers_mol.Atoms)

    pers_atoms_for_clipping: FSharpList[AtomInfo] = to_list(delay(_arrow47))
    def _arrow49(__unit: None=None, pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> IEnumerable_1[Tuple[Point3D, float, AtomInfo]]:
        def _arrow48(match_value: Tuple[AtomInfo, AtomInfo]) -> IEnumerable_1[Tuple[Point3D, float, AtomInfo]]:
            pers_other_atom: AtomInfo = match_value[0]
            if AtomInfo__Intersects_Z6467FD91(atom, match_value[1]):
                match_value_2: SphereSphereIntersection = AtomInfo__Intersection_Z6467FD91(pers_atom, pers_other_atom)
                if match_value_2.tag == 3:
                    p: Point3D = match_value_2.fields[0]
                    return singleton((p, match_value_2.fields[1] * 0.8, pers_other_atom)) if (Point3D__Distance_591E286D(pov, p) < dist_pov_atom) else empty()

                else: 
                    return empty()


            else: 
                return empty()


        return collect(_arrow48, zip(pers_atoms_for_clipping, atoms_for_clipping))

    clipping_atoms: Array[Tuple[Point3D, float, AtomInfo]] = to_array(delay(_arrow49))
    if len(clipping_atoms) == 0:
        return empty_1()

    else: 
        def mapping_1(tupled_arg: Tuple[Point3D, float, AtomInfo], pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> Tuple[Optional[Tuple[Point2D, Point2D]], Point2D, float, Point2D, float, AtomInfo, float]:
            i_p: Point3D = tupled_arg[0]
            i_r: float = tupled_arg[1]
            other: AtomInfo = tupled_arg[2]
            p1: Point2D = Point2D(pers_atom.Center.X, pers_atom.Center.Y)
            p2: Point2D = Point2D(i_p.X, i_p.Y)
            return (intersection_between_circles(p1, pers_atom.Radius, p2, i_r), p1, pers_atom.Radius, p2, i_r, other, other.Radius)

        intersections: Array[Tuple[Optional[Tuple[Point2D, Point2D]], Point2D, float, Point2D, float, AtomInfo, float]] = map_1(mapping_1, clipping_atoms, None)
        def _arrow52(__unit: None=None, pov: Point3D=pov, pers_atom: AtomInfo=pers_atom, pers_mol: Molecule=pers_mol, atom: AtomInfo=atom, mol: Molecule=mol) -> IEnumerable_1[ClipPath]:
            def _arrow51(match_value_3: Tuple[Optional[Tuple[Point2D, Point2D]], Point2D, float, Point2D, float, AtomInfo, float]) -> IEnumerable_1[ClipPath]:
                r1: float = match_value_3[2]
                p2_1: Point2D = match_value_3[3]
                p1_1: Point2D = match_value_3[1]
                match_value_4: Optional[Tuple[Point2D, Point2D]] = match_value_3[0]
                if match_value_4 is not None:
                    l2: Point2D = match_value_4[1]
                    l1: Point2D = match_value_4[0]
                    def _arrow50(__unit: None=None) -> SelectForSide:
                        match_value_6: float = Point2D__Distance_591E284C(p1_1, p2_1)
                        return SelectForSide(2) if ((r1 < match_value_3[6]) if (match_value_6 < r1) else False) else (SelectForSide(0, p1_1) if (match_value_6 < r1) else SelectForSide(1, p1_1))

                    return singleton(ClipPath((l1, l2), _arrow50() if same_side_of_line(l1, l2, p1_1, p2_1) else SelectForSide(0, p1_1)))

                else: 
                    return empty()


            return collect(_arrow51, intersections)

        return to_list(delay(_arrow52))



def intersection_between_lines(l1_: Point2D, l1__1: Point2D, l2_: Point2D, l2__1: Point2D) -> Optional[Point2D]:
    l1: Tuple[Point2D, Point2D] = (l1_, l1__1)
    l2: Tuple[Point2D, Point2D] = (l2_, l2__1)
    p1: Point2D = l1[0]
    p3: Point2D = l2[0]
    a1: float = calc_slope(p1, l1[1])
    a2: float = calc_slope(p3, l2[1])
    if a1 == a2:
        return None

    else: 
        c1: float = p1.Y - (a1 * p1.X)
        c2: float = p3.Y - (a2 * p3.X)
        if True if (True if (c1 == inf) else (c1 == (-inf))) else (True if (c2 == inf) else (c2 == (-inf))):
            return None

        else: 
            x: float = divide(c1 - c2, a2 - a1)
            return Point2D(x, (a1 * x) + c1)




__all__ = ["intersection_between_circles", "calc_slope", "same_side_of_line", "clip", "intersection_between_lines"]

