from __future__ import annotations
from collections.abc import Callable
from dataclasses import dataclass
from math import (pow, cos, sin)
from typing import Any
from .cine_mol_helpers import clamp
from .fable_modules.fable_library.double import (divide, sqrt)
from .fable_modules.fable_library.list import (length as length_1, fold, FSharpList, map as map_1, find, filter, max)
from .fable_modules.fable_library.range import range_big_int
from .fable_modules.fable_library.reflection import (TypeInfo, int32_type, union_type, float64_type, record_type, option_type, list_type, string_type, bool_type)
from .fable_modules.fable_library.seq import (to_list, delay, map, collect, empty, singleton)
from .fable_modules.fable_library.string_ import (to_text, interpolate, join)
from .fable_modules.fable_library.types import (Array, Union, Record, to_string)
from .fable_modules.fable_library.util import (IEnumerable_1, equals, compare_primitives)

def _expr48() -> TypeInfo:
    return union_type("CineMol.Types.Color", [], Color, lambda: [[("Item1", int32_type), ("Item2", int32_type), ("Item3", int32_type)]])


class Color(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Color"]

    def __str__(self, __unit: None=None) -> str:
        this: Color = self
        return ((((("rgb(" + str(this.fields[0])) + ",") + str(this.fields[1])) + ",") + str(this.fields[2])) + ")"


Color_reflection = _expr48

def Color__get_ToHex(this: Color) -> str:
    return "#%P(x2)%P(x2)%P(x2)"


def Color__Diffuse_5E38073B(this: Color, alpha: float) -> Color:
    alpha_1: float = clamp(0.0, 1.0, alpha)
    def diffuse_channel(c: int, this: Any=this, alpha: Any=alpha) -> int:
        return int(c * alpha_1)

    tupled_arg: tuple[int, int, int] = (diffuse_channel(this.fields[0]), diffuse_channel(this.fields[1]), diffuse_channel(this.fields[2]))
    return Color(0, tupled_arg[0], tupled_arg[1], tupled_arg[2])


def _expr49() -> TypeInfo:
    return union_type("CineMol.Types.ModelStyle", [], ModelStyle, lambda: [[], [], []])


class ModelStyle(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["SpaceFilling", "BallAndStick", "WireFrame"]


ModelStyle_reflection = _expr49

def _expr50() -> TypeInfo:
    return union_type("CineMol.Types.ArtStyle", [], ArtStyle, lambda: [[], []])


class ArtStyle(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Cartoon", "Glossy"]


ArtStyle_reflection = _expr50

def _expr51() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Point2D", [], Geometry_Point2D, lambda: [("X", float64_type), ("Y", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Point2D(Record):
    X: float
    Y: float

Geometry_Point2D_reflection = _expr51

def _expr52() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Vector2D", [], Geometry_Vector2D, lambda: [("X", float64_type), ("Y", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Vector2D(Record):
    X: float
    Y: float

Geometry_Vector2D_reflection = _expr52

def _expr53() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Point3D", [], Geometry_Point3D, lambda: [("X", float64_type), ("Y", float64_type), ("Z", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Point3D(Record):
    X: float
    Y: float
    Z: float

Geometry_Point3D_reflection = _expr53

def _expr54() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Vector3D", [], Geometry_Vector3D, lambda: [("X", float64_type), ("Y", float64_type), ("Z", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Vector3D(Record):
    X: float
    Y: float
    Z: float

Geometry_Vector3D_reflection = _expr54

def _expr55() -> TypeInfo:
    return union_type("CineMol.Types.Geometry.Axis", [], Geometry_Axis, lambda: [[], [], []])


class Geometry_Axis(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["X", "Y", "Z"]


Geometry_Axis_reflection = _expr55

def _expr56() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Circle", [], Geometry_Circle, lambda: [("Center", Geometry_Point2D_reflection()), ("Radius", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Circle(Record):
    Center: Geometry_Point2D
    Radius: float

Geometry_Circle_reflection = _expr56

def _expr57() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Sphere", [], Geometry_Sphere, lambda: [("Center", Geometry_Point3D_reflection()), ("Radius", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Sphere(Record):
    Center: Geometry_Point3D
    Radius: float

Geometry_Sphere_reflection = _expr57

def Geometry_Point2D_op_Addition_50440620(p1: Geometry_Point2D, p2: Geometry_Point2D) -> Geometry_Point2D:
    return Geometry_Point2D(p1.X + p2.X, p1.Y + p2.Y)


def Geometry_Point2D_op_Subtraction_50440620(p1: Geometry_Point2D, p2: Geometry_Point2D) -> Geometry_Point2D:
    return Geometry_Point2D(p1.X - p2.X, p1.Y - p2.Y)


def Geometry_Point2D_op_Multiply_50440620(p1: Geometry_Point2D, p2: Geometry_Point2D) -> Geometry_Point2D:
    return Geometry_Point2D(p1.X * p2.X, p1.Y * p2.Y)


def Geometry_Point2D__Add_5E38073B(this: Geometry_Point2D, v: float) -> Geometry_Point2D:
    return Geometry_Point2D(this.X + v, this.Y + v)


def Geometry_Point2D__Mul_5E38073B(this: Geometry_Point2D, v: float) -> Geometry_Point2D:
    return Geometry_Point2D(this.X * v, this.X * v)


def Geometry_Point2D__Div_5E38073B(this: Geometry_Point2D, v: float) -> Geometry_Point2D:
    return Geometry_Point2D(divide(this.X, v), divide(this.Y, v))


def Geometry_Point2D__Pow_5E38073B(this: Geometry_Point2D, v: float) -> Geometry_Point2D:
    return Geometry_Point2D(pow(this.X, v), pow(this.Y, v))


def Geometry_Point2D__Sum(this: Geometry_Point2D) -> float:
    return this.X + this.Y


def Geometry_Point2D__Dist_2282202F(this: Geometry_Point2D, other: Geometry_Point2D) -> float:
    return sqrt(Geometry_Point2D__Sum(Geometry_Point2D__Pow_5E38073B(Geometry_Point2D_op_Subtraction_50440620(this, other), 2.0)))


def Geometry_Point2D__Midpoint_2282202F(this: Geometry_Point2D, other: Geometry_Point2D) -> Geometry_Point2D:
    return Geometry_Point2D__Div_5E38073B(Geometry_Point2D_op_Addition_50440620(this, other), 2.0)


def Geometry_Point2D__CreateVector_2282202F(this: Geometry_Point2D, other: Geometry_Point2D) -> Geometry_Vector2D:
    return Geometry_Vector2D(other.X - this.X, other.Y - this.Y)


def Geometry_Point2D_Centroid_18E940CD(ps: FSharpList[Geometry_Point2D]) -> Geometry_Point2D:
    def folder(p_sum: Geometry_Point2D, p: Geometry_Point2D, ps: Any=ps) -> Geometry_Point2D:
        return Geometry_Point2D_op_Addition_50440620(p_sum, p)

    return Geometry_Point2D__Div_5E38073B(fold(folder, Geometry_Point2D(0.0, 0.0), ps), length_1(ps))


def Geometry_Vector2D__Cross_66BEF79A(this: Geometry_Vector2D, other: Geometry_Vector2D) -> float:
    return (this.X * other.Y) - (this.Y * other.X)


def Geometry_Point3D_op_Addition_504401C0(p1: Geometry_Point3D, p2: Geometry_Point3D) -> Geometry_Point3D:
    return Geometry_Point3D(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z)


def Geometry_Point3D_op_Subtraction_504401C0(p1: Geometry_Point3D, p2: Geometry_Point3D) -> Geometry_Point3D:
    return Geometry_Point3D(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z)


def Geometry_Point3D_op_Multiply_504401C0(p1: Geometry_Point3D, p2: Geometry_Point3D) -> Geometry_Point3D:
    return Geometry_Point3D(p1.X * p2.X, p1.Y * p2.Y, p1.Z * p2.Z)


def Geometry_Point3D__Add_5E38073B(this: Geometry_Point3D, v: float) -> Geometry_Point3D:
    return Geometry_Point3D(this.X + v, this.Y + v, this.Z + v)


def Geometry_Point3D__Mul_5E38073B(this: Geometry_Point3D, v: float) -> Geometry_Point3D:
    return Geometry_Point3D(this.X * v, this.Y * v, this.Z * v)


def Geometry_Point3D__Div_5E38073B(this: Geometry_Point3D, v: float) -> Geometry_Point3D:
    return Geometry_Point3D(divide(this.X, v), divide(this.Y, v), divide(this.Z, v))


def Geometry_Point3D__Pow_5E38073B(this: Geometry_Point3D, v: float) -> Geometry_Point3D:
    return Geometry_Point3D(pow(this.X, v), pow(this.Y, v), pow(this.Z, v))


def Geometry_Point3D__Sum(this: Geometry_Point3D) -> float:
    return (this.X + this.Y) + this.Z


def Geometry_Point3D__Dist_2282200E(this: Geometry_Point3D, other: Geometry_Point3D) -> float:
    return sqrt(Geometry_Point3D__Sum(Geometry_Point3D__Pow_5E38073B(Geometry_Point3D_op_Subtraction_504401C0(this, other), 2.0)))


def Geometry_Point3D__Midpoint_2282200E(this: Geometry_Point3D, other: Geometry_Point3D) -> Geometry_Point3D:
    return Geometry_Point3D__Div_5E38073B(Geometry_Point3D_op_Addition_504401C0(this, other), 2.0)


def Geometry_Point3D__VectorTo_2282200E(this: Geometry_Point3D, other: Geometry_Point3D) -> Geometry_Vector3D:
    return Geometry_Vector3D(other.X - this.X, other.Y - this.Y, other.Z - this.Z)


def Geometry_Point3D__MoveTowards(this: Geometry_Point3D, other: Geometry_Point3D, distance: float) -> Geometry_Point3D:
    match_value: float = Geometry_Point3D__Dist_2282200E(this, other) - distance
    if match_value <= 0.0:
        return other

    else: 
        new_dist_2: float = match_value
        norm: Geometry_Point3D = Geometry_Vector3D__Normalize(Geometry_Point3D__VectorTo_2282200E(this, other))
        return Geometry_Point3D_op_Addition_504401C0(this, Geometry_Point3D(norm.X * new_dist_2, norm.Y * new_dist_2, norm.Z * new_dist_2))



def Geometry_Point3D__Rotate(p: Geometry_Point3D, axis: Geometry_Axis, rad: float) -> Geometry_Point3D:
    if axis.tag == 1:
        return Geometry_Point3D((p.X * cos(rad)) + (p.Z * sin(rad)), p.Y, (p.Z * cos(rad)) - (p.X * sin(rad)))

    elif axis.tag == 2:
        return Geometry_Point3D((p.X * cos(rad)) - (p.Y * sin(rad)), (p.X * sin(rad)) + (p.Y * cos(rad)), p.Z)

    else: 
        return Geometry_Point3D(p.X, (p.Y * cos(rad)) - (p.Z * sin(rad)), (p.Y * sin(rad)) + (p.Z * cos(rad)))



def Geometry_Point3D__ToPoint2D(this: Geometry_Point3D) -> Geometry_Vector2D:
    return Geometry_Vector2D(this.X, this.Y)


def Geometry_Point3D_Centroid_Z3431E334(ps: FSharpList[Geometry_Point3D]) -> Geometry_Point3D:
    def folder(p_sum: Geometry_Point3D, p: Geometry_Point3D, ps: Any=ps) -> Geometry_Point3D:
        return Geometry_Point3D_op_Addition_504401C0(p_sum, p)

    return Geometry_Point3D__Div_5E38073B(fold(folder, Geometry_Point3D(0.0, 0.0, 0.0), ps), length_1(ps))


def Geometry_Vector3D__Normalize(this: Geometry_Vector3D) -> Geometry_Point3D:
    length: float = sqrt(((this.X * this.X) + (this.Y * this.Y)) + (this.Z * this.Z))
    return Geometry_Point3D(divide(this.X, length), divide(this.Y, length), divide(this.Z, length))


def Geometry_Axis_Origin(__unit: None=None) -> Geometry_Point3D:
    return Geometry_Point3D(0.0, 0.0, 0.0)


def Geometry_Sphere__Encloses_Z460483F4(this: Geometry_Sphere, o: Geometry_Sphere) -> bool:
    return (Geometry_Point3D__Dist_2282200E(this.Center, o.Center) + o.Radius) <= this.Radius


def Geometry_Sphere__Intersects_Z460483F4(this: Geometry_Sphere, o: Geometry_Sphere) -> bool:
    return Geometry_Point3D__Dist_2282200E(this.Center, o.Center) <= (this.Radius + o.Radius)


def Geometry_Sphere__PointsOnSphere_Z524259A4(this: Geometry_Sphere, resolution: int) -> FSharpList[Geometry_Point3D]:
    N: float = resolution
    num_points_phi: int = int(divide(N, 2.0)) or 0
    def _arrow59(__unit: None=None, this: Any=this, resolution: Any=resolution) -> IEnumerable_1[float]:
        def _arrow58(i: int) -> float:
            return divide(i, num_points_phi - 1.0) * 3.141592653589793

        return map(_arrow58, to_list(range_big_int(1, 1, num_points_phi)))

    phis: FSharpList[float] = to_list(delay(_arrow59))
    num_points_theta: int = int(divide(N, 2.0)) or 0
    def _arrow61(__unit: None=None, this: Any=this, resolution: Any=resolution) -> IEnumerable_1[float]:
        def _arrow60(i_1: int) -> float:
            return (divide(i_1, num_points_theta) * 2.0) * 3.141592653589793

        return map(_arrow60, to_list(range_big_int(0, 1, num_points_theta)))

    thetas: FSharpList[float] = to_list(delay(_arrow61))
    def _arrow64(__unit: None=None, this: Any=this, resolution: Any=resolution) -> IEnumerable_1[Geometry_Point3D]:
        def _arrow63(phi: float) -> IEnumerable_1[Geometry_Point3D]:
            z: float = this.Center.Z + (this.Radius * cos(phi))
            if z < this.Center.Z:
                return empty()

            else: 
                def _arrow62(theta: float) -> IEnumerable_1[Geometry_Point3D]:
                    return singleton(Geometry_Point3D(this.Center.X + ((this.Radius * sin(phi)) * cos(theta)), this.Center.Y + ((this.Radius * sin(phi)) * sin(theta)), z))

                return collect(_arrow62, thetas)


        return collect(_arrow63, phis)

    return to_list(delay(_arrow64))


def _expr65() -> TypeInfo:
    return record_type("CineMol.Types.Geometry.Cylinder", [], Geometry_Cylinder, lambda: [("Start", Geometry_Point3D_reflection()), ("End", Geometry_Point3D_reflection()), ("Radius", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Geometry_Cylinder(Record):
    Start: Geometry_Point3D
    End: Geometry_Point3D
    Radius: float

Geometry_Cylinder_reflection = _expr65

def Geometry_Cylinder__IsInside_2282200E(this: Geometry_Cylinder, point: Geometry_Point3D) -> bool:
    cylinder_direction: Geometry_Point3D = Geometry_Point3D_op_Subtraction_504401C0(this.End, this.Start)
    def _arrow66(__unit: None=None, this: Any=this, point: Any=point) -> float:
        p: Geometry_Point3D = Geometry_Point3D_op_Subtraction_504401C0(point, this.Start)
        return ((p.X * cylinder_direction.X) + (p.Y * cylinder_direction.Y)) + (p.Z * cylinder_direction.Z)

    def _arrow67(__unit: None=None, this: Any=this, point: Any=point) -> float:
        c: Geometry_Point3D = cylinder_direction
        return ((c.X * c.X) + (c.Y * c.Y)) + (c.Z * c.Z)

    projection: float = divide(_arrow66(), _arrow67())
    if projection < 0.0:
        return False

    elif projection > 1.0:
        return False

    else: 
        def _arrow68(__unit: None=None, this: Any=this, point: Any=point) -> float:
            d: Geometry_Point3D = Geometry_Point3D_op_Subtraction_504401C0(point, Geometry_Point3D_op_Multiply_504401C0(Geometry_Point3D__Add_5E38073B(this.Start, projection), cylinder_direction))
            return ((d.X * d.X) + (d.Y * d.Y)) + (d.Z * d.Z)

        return _arrow68() <= (this.Radius * this.Radius)



def Geometry_Cylinder__PointsOnCylinder_Z524259A4(this: Geometry_Cylinder, resolution: int) -> FSharpList[Geometry_Point3D]:
    N: int = resolution or 0
    latitude_divisions: int = (N // 2) or 0
    longitude_divisions: int = N or 0
    def _arrow71(__unit: None=None, this: Any=this, resolution: Any=resolution) -> IEnumerable_1[Geometry_Point3D]:
        def _arrow70(lat: int) -> IEnumerable_1[Geometry_Point3D]:
            theta: float = divide(lat * 3.141592653589793, latitude_divisions)
            sin_theta: float = sin(theta)
            cos_theta: float = cos(theta)
            def _arrow69(lon: int) -> IEnumerable_1[Geometry_Point3D]:
                phi: float = divide((lon * 2.0) * 3.141592653589793, longitude_divisions)
                sin_phi: float = sin(phi)
                return singleton(Geometry_Point3D((this.Radius * sin_theta) * cos(phi), (this.Radius * sin_theta) * sin_phi, this.Radius * cos_theta))

            return collect(_arrow69, range_big_int(0, 1, longitude_divisions))

        return collect(_arrow70, range_big_int(0, 1, latitude_divisions))

    return to_list(delay(_arrow71))


def Geometry_calcCentroid(points: FSharpList[Geometry_Point2D]) -> Geometry_Point2D:
    def folder(p_sum: Geometry_Point2D, p: Geometry_Point2D, points: Any=points) -> Geometry_Point2D:
        return Geometry_Point2D_op_Addition_50440620(p_sum, p)

    return Geometry_Point2D__Div_5E38073B(fold(folder, Geometry_Point2D(0.0, 0.0), points), length_1(points))


def _expr72() -> TypeInfo:
    return union_type("CineMol.Types.Chem.AtomType", [], Chem_AtomType, lambda: [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []])


class Chem_AtomType(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra"]


Chem_AtomType_reflection = _expr72

def Chem_AtomType_FromString_Z721C83C5(atom_string: str) -> Chem_AtomType | None:
    if atom_string == "H":
        return Chem_AtomType(0)

    elif atom_string == "He":
        return Chem_AtomType(1)

    elif atom_string == "Li":
        return Chem_AtomType(2)

    elif atom_string == "Be":
        return Chem_AtomType(3)

    elif atom_string == "B":
        return Chem_AtomType(4)

    elif atom_string == "C":
        return Chem_AtomType(5)

    elif atom_string == "N":
        return Chem_AtomType(6)

    elif atom_string == "O":
        return Chem_AtomType(7)

    elif atom_string == "F":
        return Chem_AtomType(8)

    elif atom_string == "Ne":
        return Chem_AtomType(9)

    elif atom_string == "Na":
        return Chem_AtomType(10)

    elif atom_string == "Mg":
        return Chem_AtomType(11)

    elif atom_string == "Al":
        return Chem_AtomType(12)

    elif atom_string == "Si":
        return Chem_AtomType(13)

    elif atom_string == "P":
        return Chem_AtomType(14)

    elif atom_string == "S":
        return Chem_AtomType(15)

    elif atom_string == "Cl":
        return Chem_AtomType(16)

    elif atom_string == "Ar":
        return Chem_AtomType(17)

    elif atom_string == "K":
        return Chem_AtomType(18)

    elif atom_string == "Ca":
        return Chem_AtomType(19)

    elif atom_string == "Sc":
        return Chem_AtomType(20)

    elif atom_string == "Ti":
        return Chem_AtomType(21)

    elif atom_string == "V":
        return Chem_AtomType(22)

    elif atom_string == "Cr":
        return Chem_AtomType(23)

    elif atom_string == "Mn":
        return Chem_AtomType(24)

    elif atom_string == "Fe":
        return Chem_AtomType(25)

    elif atom_string == "Co":
        return Chem_AtomType(26)

    elif atom_string == "Ni":
        return Chem_AtomType(27)

    elif atom_string == "Cu":
        return Chem_AtomType(28)

    elif atom_string == "Zn":
        return Chem_AtomType(29)

    elif atom_string == "Ga":
        return Chem_AtomType(30)

    elif atom_string == "Ge":
        return Chem_AtomType(31)

    elif atom_string == "As":
        return Chem_AtomType(32)

    elif atom_string == "Se":
        return Chem_AtomType(33)

    elif atom_string == "Br":
        return Chem_AtomType(34)

    elif atom_string == "Kr":
        return Chem_AtomType(35)

    elif atom_string == "Rb":
        return Chem_AtomType(36)

    elif atom_string == "Sr":
        return Chem_AtomType(37)

    elif atom_string == "Y":
        return Chem_AtomType(38)

    elif atom_string == "Zr":
        return Chem_AtomType(39)

    elif atom_string == "Nb":
        return Chem_AtomType(40)

    elif atom_string == "Mo":
        return Chem_AtomType(41)

    elif atom_string == "Tc":
        return Chem_AtomType(42)

    elif atom_string == "Ru":
        return Chem_AtomType(43)

    elif atom_string == "Rh":
        return Chem_AtomType(44)

    elif atom_string == "Pd":
        return Chem_AtomType(45)

    elif atom_string == "Ag":
        return Chem_AtomType(46)

    elif atom_string == "Cd":
        return Chem_AtomType(47)

    elif atom_string == "In":
        return Chem_AtomType(48)

    elif atom_string == "Sn":
        return Chem_AtomType(49)

    elif atom_string == "Sb":
        return Chem_AtomType(50)

    elif atom_string == "Te":
        return Chem_AtomType(51)

    elif atom_string == "I":
        return Chem_AtomType(52)

    elif atom_string == "Xe":
        return Chem_AtomType(53)

    elif atom_string == "Cs":
        return Chem_AtomType(54)

    elif atom_string == "Ba":
        return Chem_AtomType(55)

    elif atom_string == "Lu":
        return Chem_AtomType(56)

    elif atom_string == "Hf":
        return Chem_AtomType(57)

    elif atom_string == "Ta":
        return Chem_AtomType(58)

    elif atom_string == "W":
        return Chem_AtomType(59)

    elif atom_string == "Re":
        return Chem_AtomType(60)

    elif atom_string == "Os":
        return Chem_AtomType(61)

    elif atom_string == "Ir":
        return Chem_AtomType(62)

    elif atom_string == "Pt":
        return Chem_AtomType(63)

    elif atom_string == "Au":
        return Chem_AtomType(64)

    elif atom_string == "Hg":
        return Chem_AtomType(65)

    elif atom_string == "Tl":
        return Chem_AtomType(66)

    elif atom_string == "Pb":
        return Chem_AtomType(67)

    elif atom_string == "Bi":
        return Chem_AtomType(68)

    elif atom_string == "Po":
        return Chem_AtomType(69)

    elif atom_string == "At":
        return Chem_AtomType(70)

    elif atom_string == "Rn":
        return Chem_AtomType(71)

    elif atom_string == "Fr":
        return Chem_AtomType(72)

    elif atom_string == "Ra":
        return Chem_AtomType(73)

    else: 
        return None



def _expr73() -> TypeInfo:
    return union_type("CineMol.Types.Chem.BondType", [], Chem_BondType, lambda: [[], [], [], []])


class Chem_BondType(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Single", "Double", "Triple", "Aromatic"]


Chem_BondType_reflection = _expr73

def Chem_BondType_FromString_Z721C83C5(bond_string: str) -> Chem_BondType | None:
    (pattern_matching_result,) = (None,)
    if bond_string == "1":
        pattern_matching_result = 0

    elif bond_string == "SINGLE":
        pattern_matching_result = 0

    elif bond_string == "Single":
        pattern_matching_result = 0

    elif bond_string == "2":
        pattern_matching_result = 1

    elif bond_string == "DOUBLE":
        pattern_matching_result = 1

    elif bond_string == "Double":
        pattern_matching_result = 1

    elif bond_string == "3":
        pattern_matching_result = 2

    elif bond_string == "TRIPLE":
        pattern_matching_result = 2

    elif bond_string == "Triple":
        pattern_matching_result = 2

    elif bond_string == "4":
        pattern_matching_result = 3

    elif bond_string == "AROMATIC":
        pattern_matching_result = 3

    elif bond_string == "Aromatic":
        pattern_matching_result = 3

    else: 
        pattern_matching_result = 4

    if pattern_matching_result == 0:
        return Chem_BondType(0)

    elif pattern_matching_result == 1:
        return Chem_BondType(1)

    elif pattern_matching_result == 2:
        return Chem_BondType(2)

    elif pattern_matching_result == 3:
        return Chem_BondType(3)

    elif pattern_matching_result == 4:
        return None



def _expr74() -> TypeInfo:
    return record_type("CineMol.Types.Chem.Atom", [], Chem_Atom, lambda: [("Index", int32_type), ("Type", Chem_AtomType_reflection()), ("Color", Color_reflection()), ("Opacity", float64_type), ("Position", Geometry_Point3D_reflection()), ("Radius", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Chem_Atom(Record):
    Index: int
    Type: Chem_AtomType
    Color: Color
    Opacity: float
    Position: Geometry_Point3D
    Radius: float

Chem_Atom_reflection = _expr74

def _expr75() -> TypeInfo:
    return record_type("CineMol.Types.Chem.Bond", [], Chem_Bond, lambda: [("BeginIndex", int32_type), ("EndIndex", int32_type), ("Type", Chem_BondType_reflection()), ("BeginAtomIndex", int32_type), ("EndAtomIndex", int32_type), ("Opacity", option_type(float64_type)), ("Color", option_type(Color_reflection())), ("Radius", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Chem_Bond(Record):
    BeginIndex: int
    EndIndex: int
    Type: Chem_BondType
    BeginAtomIndex: int
    EndAtomIndex: int
    Opacity: float | None
    Color: Color | None
    Radius: float

Chem_Bond_reflection = _expr75

def _expr76() -> TypeInfo:
    return record_type("CineMol.Types.Chem.Molecule", [], Chem_Molecule, lambda: [("Atoms", list_type(Chem_Atom_reflection())), ("Bonds", list_type(Chem_Bond_reflection()))])


@dataclass(eq = False, repr = False, slots = True)
class Chem_Molecule(Record):
    Atoms: FSharpList[Chem_Atom]
    Bonds: FSharpList[Chem_Bond]

Chem_Molecule_reflection = _expr76

def Chem_Molecule__AdjustForCentroid(this: Chem_Molecule) -> Chem_Molecule:
    def mapping(atom: Chem_Atom, this: Any=this) -> Geometry_Point3D:
        return atom.Position

    centroid: Geometry_Point3D = Geometry_Point3D_Centroid_Z3431E334(map_1(mapping, this.Atoms))
    def mapping_1(atom_1: Chem_Atom, this: Any=this) -> Chem_Atom:
        return Chem_Atom(atom_1.Index, atom_1.Type, atom_1.Color, atom_1.Opacity, Geometry_Point3D_op_Subtraction_504401C0(atom_1.Position, centroid), atom_1.Radius)

    return Chem_Molecule(map_1(mapping_1, this.Atoms), this.Bonds)


def Chem_Molecule__GetAtom_Z524259A4(this: Chem_Molecule, atom_index: int) -> Chem_Atom | None:
    try: 
        def predicate(atom: Chem_Atom) -> bool:
            return atom.Index == atom_index

        return find(predicate, this.Atoms)

    except Exception as match_value:
        return None



def Chem_Molecule__GetBonds(this: Chem_Molecule, include_hydrogen_atoms: bool, atom_index: int) -> FSharpList[Chem_Bond]:
    def predicate_1(bond_1: Chem_Bond, this: Any=this, include_hydrogen_atoms: Any=include_hydrogen_atoms, atom_index: Any=atom_index) -> bool:
        begin_atom: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(this, bond_1.BeginAtomIndex)
        end_atom: Chem_Atom | None = Chem_Molecule__GetAtom_Z524259A4(this, bond_1.EndAtomIndex)
        (pattern_matching_result, begin_atom_1, end_atom_1) = (None, None, None)
        if begin_atom is not None:
            if end_atom is not None:
                pattern_matching_result = 0
                begin_atom_1 = begin_atom
                end_atom_1 = end_atom

            else: 
                pattern_matching_result = 1


        else: 
            pattern_matching_result = 1

        if pattern_matching_result == 0:
            if include_hydrogen_atoms:
                return True

            elif not equals(begin_atom_1.Type, Chem_AtomType(0)):
                return not equals(end_atom_1.Type, Chem_AtomType(0))

            else: 
                return False


        elif pattern_matching_result == 1:
            return False


    def predicate(bond: Chem_Bond, this: Any=this, include_hydrogen_atoms: Any=include_hydrogen_atoms, atom_index: Any=atom_index) -> bool:
        if bond.BeginAtomIndex == atom_index:
            return True

        else: 
            return bond.EndAtomIndex == atom_index


    return filter(predicate_1, filter(predicate, this.Bonds))


def _expr77() -> TypeInfo:
    return record_type("CineMol.Types.Svg.ViewBox", [], Svg_ViewBox, lambda: [("MinX", float64_type), ("MinY", float64_type), ("Width", float64_type), ("Height", float64_type)])


@dataclass(eq = False, repr = False, slots = True)
class Svg_ViewBox(Record):
    MinX: float
    MinY: float
    Width: float
    Height: float
    def __str__(self, __unit: None=None) -> str:
        this: Svg_ViewBox = self
        return to_text(interpolate("viewBox=\"%.3f%P() %.3f%P() %.3f%P() %.3f%P()\"", [this.MinX, this.MinY, this.Width, this.Height]))


Svg_ViewBox_reflection = _expr77

def _expr78() -> TypeInfo:
    return union_type("CineMol.Types.Svg.Fill", [], Svg_Fill, lambda: [[("Item1", int32_type), ("Item2", Geometry_Point2D_reflection()), ("Item3", float64_type), ("Item4", Color_reflection())], [("Item1", int32_type), ("Item2", Geometry_Point2D_reflection()), ("Item3", Geometry_Point2D_reflection()), ("Item4", Color_reflection())]])


class Svg_Fill(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["RadialGradient", "LinearGradient"]

    def __str__(self, __unit: None=None) -> str:
        this: Svg_Fill = self
        if this.tag == 0:
            color_1: Color = this.fields[3]
            center: Geometry_Point2D = this.fields[1]
            cx: float = center.X
            cy: float = center.Y
            return to_text(interpolate("<radialGradient id=\"item-%P()\" cx=\"%.3f%P()\" cy=\"%.3f%P()\" r=\"%.3f%P()\" fx=\"%.3f%P()\" fy=\"%.3f%P()\" gradientTransform=\"matrix(1,0,0,1,0,0)\" gradientUnits=\"userSpaceOnUse\"><stop offset=\"0.00\" stop-color=\"%P()\"/><stop offset=\"1.00\" stop-color=\"%P()\"/></radialGradient>", [this.fields[0], cx, cy, this.fields[2] * 1.5, cx, cy, color_1, Color__Diffuse_5E38073B(color_1, 0.5)]))

        else: 
            stop: Geometry_Point2D = this.fields[2]
            start: Geometry_Point2D = this.fields[1]
            color: Color = this.fields[3]
            return to_text(interpolate("<linearGradient id=\"item-%P()\" x1=\"%.3f%P()\" x2=\"%.3f%P()\" y1=\"%.3f%P()\" y2=\"%.3f%P()\" gradientUnits=\"userSpaceOnUse\" spreadMethod=\"reflect\"><stop offset=\"0.00\" stop-color=\"%P()\"/><stop offset=\"1.00\" stop-color=\"%P()\"/></linearGradient>", [this.fields[0], start.X, stop.X, start.Y, stop.Y, color, Color__Diffuse_5E38073B(color, 0.5)]))



Svg_Fill_reflection = _expr78

def _expr79() -> TypeInfo:
    return union_type("CineMol.Types.Svg.Shape", [], Svg_Shape, lambda: [[("Item1", int32_type), ("Item2", Color_reflection()), ("Item3", Geometry_Circle_reflection()), ("Item4", float64_type)], [("Item1", int32_type), ("Item2", Color_reflection()), ("Item3", list_type(Geometry_Point2D_reflection())), ("Item4", float64_type)], [("Item1", int32_type), ("Item2", Color_reflection()), ("Item3", Geometry_Point2D_reflection()), ("Item4", Geometry_Point2D_reflection()), ("Item5", Geometry_Point2D_reflection()), ("Item6", Geometry_Point2D_reflection()), ("Item7", float64_type)], [("Item1", int32_type), ("Item2", Color_reflection()), ("Item3", Geometry_Point2D_reflection()), ("Item4", Geometry_Point2D_reflection()), ("Item5", float64_type)]])


class Svg_Shape(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Circle", "Polygon", "BondPolygon", "Line"]


Svg_Shape_reflection = _expr79

def Svg_Shape__ToSvg_DEA3CF1(this: Svg_Shape, style: ArtStyle) -> str:
    if this.tag == 1:
        opacity_1: float = this.fields[3]
        def mapping(p: Geometry_Point2D, this: Any=this, style: Any=style) -> str:
            return to_text(interpolate("%.3f%P(),%.3f%P()", [p.X, p.Y]))

        points_str: str = join(" ", map_1(mapping, this.fields[2]))
        if style.tag == 1:
            return ((((("<polygon class=\"item-" + str(this.fields[0])) + "\" points=\"") + points_str) + "\" fill-opacity=\"") + str(opacity_1)) + "\"/>"

        else: 
            return ((((((("<polygon points=\"" + points_str) + "\" fill=\"") + str(this.fields[1])) + "\" style=\"stroke:black;stroke-width:0.05\" stroke-linejoin=\"round\" fill-opacity=\"") + str(opacity_1)) + "\" stroke-opacity=\"") + str(opacity_1)) + "\"/>"


    elif this.tag == 2:
        p4: Geometry_Point2D = this.fields[5]
        p3: Geometry_Point2D = this.fields[4]
        p2: Geometry_Point2D = this.fields[3]
        p1: Geometry_Point2D = this.fields[2]
        opacity_2: float = this.fields[6]
        r1: float = divide(Geometry_Point2D__Dist_2282202F(p1, p2), 2.0)
        r2: float = divide(Geometry_Point2D__Dist_2282202F(p3, p4), 2.0)
        path: str = to_text(interpolate("M %.3f%P() %.3f%P() A %.3f%P() %.3f%P() 0 0 1 %.3f%P() %.3f%P() L %.3f%P() %.3f%P() A %.3f%P() %.3f%P() 0 0 1 %.3f%P() %.3f%P() Z", [p1.X, p1.Y, r1, r1, p2.X, p2.Y, p3.X, p3.Y, r2, r2, p4.X, p4.Y]))
        if style.tag == 1:
            return ((((("<path class=\"item-" + str(this.fields[0])) + "\" d=\"") + path) + "\" fill-opacity=\"") + str(opacity_2)) + "\"/>"

        else: 
            return ((((((("<path d=\"" + path) + "\" fill=\"") + str(this.fields[1])) + "\" style=\"stroke:black;stroke-width:0.05\" stroke-linejoin=\"round\" fill-opacity=\"") + str(opacity_2)) + "\" stroke-opacity=\"") + str(opacity_2)) + "\"/>"


    elif this.tag == 3:
        p2_1: Geometry_Point2D = this.fields[3]
        p1_1: Geometry_Point2D = this.fields[2]
        return to_text(interpolate("<line class=\"item-%P()\" x1=\"%.3f%P()\" y1=\"%.3f%P()\" x2=\"%.3f%P()\" y2=\"%.3f%P()\" stroke=\"%P()\" stroke-width=\"0.1\" stroke-linecap=\"round\" stroke-opacity=\"%P()\"/>", [this.fields[0], p1_1.X, p1_1.Y, p2_1.X, p2_1.Y, this.fields[1], this.fields[4]]))

    else: 
        opacity: float = this.fields[3]
        circle: Geometry_Circle = this.fields[2]
        x: float = circle.Center.X
        y: float = circle.Center.Y
        if style.tag == 1:
            return to_text(interpolate("<circle class=\"item-%P()\" cx=\"%.3f%P()\" cy=\"%.3f%P()\" r=\"%.3f%P()\" fill-opacity=\"%P()\"/>", [this.fields[0], x, y, circle.Radius, opacity]))

        else: 
            return to_text(interpolate("<circle cx=\"%.3f%P()\" cy=\"%.3f%P()\" r=\"%.3f%P()\" fill=\"%P()\" style=\"stroke:black;stroke-width:0.05\" fill-opacity=\"%P()\" stroke-opacity=\"%P()\"/>", [x, y, circle.Radius, this.fields[1], opacity, opacity]))




def Svg_Shape__Fill_DEA3CF1(this: Svg_Shape, style: ArtStyle) -> Svg_Fill | None:
    if style.tag == 1:
        if this.tag == 0:
            circle: Geometry_Circle = this.fields[2]
            def _arrow80(__unit: None=None, this: Any=this, style: Any=style) -> Svg_Fill:
                tupled_arg: tuple[int, Geometry_Point2D, float, Color] = (this.fields[0], circle.Center, circle.Radius, this.fields[1])
                return Svg_Fill(0, tupled_arg[0], tupled_arg[1], tupled_arg[2], tupled_arg[3])

            return _arrow80()

        elif this.tag == 1:
            points: FSharpList[Geometry_Point2D] = this.fields[2]
            centroid: Geometry_Point2D = Geometry_calcCentroid(points)
            def _arrow82(__unit: None=None, this: Any=this, style: Any=style) -> Svg_Fill:
                def mapping(p: Geometry_Point2D) -> float:
                    return Geometry_Point2D__Dist_2282202F(p, centroid)

                class ObjectExpr81:
                    @property
                    def Compare(self) -> Callable[[float, float], int]:
                        return compare_primitives

                tupled_arg_1: tuple[int, Geometry_Point2D, float, Color] = (this.fields[0], centroid, max(map_1(mapping, points), ObjectExpr81()), this.fields[1])
                return Svg_Fill(0, tupled_arg_1[0], tupled_arg_1[1], tupled_arg_1[2], tupled_arg_1[3])

            return _arrow82()

        elif this.tag == 2:
            def _arrow83(__unit: None=None, this: Any=this, style: Any=style) -> Svg_Fill:
                tupled_arg_2: tuple[int, Geometry_Point2D, Geometry_Point2D, Color] = (this.fields[0], Geometry_Point2D__Midpoint_2282202F(this.fields[2], this.fields[4]), Geometry_Point2D__Midpoint_2282202F(this.fields[3], this.fields[5]), this.fields[1])
                return Svg_Fill(1, tupled_arg_2[0], tupled_arg_2[1], tupled_arg_2[2], tupled_arg_2[3])

            return _arrow83()

        else: 
            return None


    else: 
        return None



def _expr84() -> TypeInfo:
    return union_type("CineMol.Types.Svg.Header", [], Svg_Header, lambda: [[("version", float64_type), ("encoding", string_type)]])


class Svg_Header(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Header"]

    def __str__(self, __unit: None=None) -> str:
        this: Svg_Header = self
        return to_text(interpolate("<?xml version=\"%.1f%P()\" encoding=\"%P()\"?>", [this.fields[0], this.fields[1]]))


Svg_Header_reflection = _expr84

def Svg_Header_New(__unit: None=None) -> Svg_Header:
    return Svg_Header(0, 1.0, "UTF-8")


def _expr85() -> TypeInfo:
    return record_type("CineMol.Types.Svg.SVG", [], Svg_SVG, lambda: [("Header", Svg_Header_reflection()), ("ID", string_type), ("ViewBox", Svg_ViewBox_reflection()), ("Objects", list_type(Svg_Shape_reflection())), ("Style", ArtStyle_reflection())])


@dataclass(eq = False, repr = False, slots = True)
class Svg_SVG(Record):
    Header: Svg_Header
    ID: str
    ViewBox: Svg_ViewBox
    Objects: FSharpList[Svg_Shape]
    Style: ArtStyle
    def __str__(self, __unit: None=None) -> str:
        this: Svg_SVG = self
        return to_string(this.Header) + Svg_SVG__Body(this)


Svg_SVG_reflection = _expr85

def Svg_SVG__Body(this: Svg_SVG) -> str:
    def mapping(x: Svg_Shape, this: Any=this) -> str:
        o: Svg_Shape = x
        (pattern_matching_result, index) = (None, None)
        if o.tag == 1:
            pattern_matching_result = 0
            index = o.fields[0]

        elif o.tag == 2:
            pattern_matching_result = 0
            index = o.fields[0]

        elif o.tag == 3:
            pattern_matching_result = 0
            index = o.fields[0]

        else: 
            pattern_matching_result = 0
            index = o.fields[0]

        if pattern_matching_result == 0:
            return (((".item-" + str(index)) + "{fill:url(#item-") + str(index)) + ");}"


    styles: str = join("\n", map_1(mapping, this.Objects))
    def mapping_1(x_1: Svg_Shape, this: Any=this) -> str:
        return Svg_Shape__ToSvg_DEA3CF1(x_1, this.Style)

    objs: str = join("\n", map_1(mapping_1, this.Objects))
    if this.Style.tag == 1:
        def mapping_2(o_1: Svg_Shape, this: Any=this) -> str:
            return to_string(Svg_Shape__Fill_DEA3CF1(o_1, this.Style))

        defs: str = join("\n", map_1(mapping_2, this.Objects))
        return ((((((((("\n<svg id=\"" + this.ID) + "\" xmlns=\"http://www.w3.org/2000/svg\" ") + to_string(this.ViewBox)) + ">\n<defs>\n<style>\n") + styles) + "\n</style>\n") + defs) + "\n</defs>\n") + objs) + "\n</svg>"

    else: 
        return ((((("<svg id=\"" + this.ID) + "\" xmlns=\"http://www.w3.org/2000/svg\" ") + to_string(this.ViewBox)) + ">\n") + objs) + "\n</svg>"



def _expr86() -> TypeInfo:
    return record_type("CineMol.Types.Drawing.DrawingOptions", [], Drawing_DrawingOptions, lambda: [("ViewBox", option_type(Svg_ViewBox_reflection())), ("ArtStyle", ArtStyle_reflection()), ("ModelStyle", ModelStyle_reflection()), ("DisplayHydrogenAtoms", bool_type), ("Resolution", int32_type)])


@dataclass(eq = False, repr = False, slots = True)
class Drawing_DrawingOptions(Record):
    ViewBox: Svg_ViewBox | None
    ArtStyle: ArtStyle
    ModelStyle: ModelStyle
    DisplayHydrogenAtoms: bool
    Resolution: int

Drawing_DrawingOptions_reflection = _expr86

def Drawing_DrawingOptions_New(__unit: None=None) -> Drawing_DrawingOptions:
    return Drawing_DrawingOptions(None, ArtStyle(0), ModelStyle(0), False, 40)


__all__ = ["Color_reflection", "Color__get_ToHex", "Color__Diffuse_5E38073B", "ModelStyle_reflection", "ArtStyle_reflection", "Geometry_Point2D_reflection", "Geometry_Vector2D_reflection", "Geometry_Point3D_reflection", "Geometry_Vector3D_reflection", "Geometry_Axis_reflection", "Geometry_Circle_reflection", "Geometry_Sphere_reflection", "Geometry_Point2D_op_Addition_50440620", "Geometry_Point2D_op_Subtraction_50440620", "Geometry_Point2D_op_Multiply_50440620", "Geometry_Point2D__Add_5E38073B", "Geometry_Point2D__Mul_5E38073B", "Geometry_Point2D__Div_5E38073B", "Geometry_Point2D__Pow_5E38073B", "Geometry_Point2D__Sum", "Geometry_Point2D__Dist_2282202F", "Geometry_Point2D__Midpoint_2282202F", "Geometry_Point2D__CreateVector_2282202F", "Geometry_Point2D_Centroid_18E940CD", "Geometry_Vector2D__Cross_66BEF79A", "Geometry_Point3D_op_Addition_504401C0", "Geometry_Point3D_op_Subtraction_504401C0", "Geometry_Point3D_op_Multiply_504401C0", "Geometry_Point3D__Add_5E38073B", "Geometry_Point3D__Mul_5E38073B", "Geometry_Point3D__Div_5E38073B", "Geometry_Point3D__Pow_5E38073B", "Geometry_Point3D__Sum", "Geometry_Point3D__Dist_2282200E", "Geometry_Point3D__Midpoint_2282200E", "Geometry_Point3D__VectorTo_2282200E", "Geometry_Point3D__MoveTowards", "Geometry_Point3D__Rotate", "Geometry_Point3D__ToPoint2D", "Geometry_Point3D_Centroid_Z3431E334", "Geometry_Vector3D__Normalize", "Geometry_Axis_Origin", "Geometry_Sphere__Encloses_Z460483F4", "Geometry_Sphere__Intersects_Z460483F4", "Geometry_Sphere__PointsOnSphere_Z524259A4", "Geometry_Cylinder_reflection", "Geometry_Cylinder__IsInside_2282200E", "Geometry_Cylinder__PointsOnCylinder_Z524259A4", "Geometry_calcCentroid", "Chem_AtomType_reflection", "Chem_AtomType_FromString_Z721C83C5", "Chem_BondType_reflection", "Chem_BondType_FromString_Z721C83C5", "Chem_Atom_reflection", "Chem_Bond_reflection", "Chem_Molecule_reflection", "Chem_Molecule__AdjustForCentroid", "Chem_Molecule__GetAtom_Z524259A4", "Chem_Molecule__GetBonds", "Svg_ViewBox_reflection", "Svg_Fill_reflection", "Svg_Shape_reflection", "Svg_Shape__ToSvg_DEA3CF1", "Svg_Shape__Fill_DEA3CF1", "Svg_Header_reflection", "Svg_Header_New", "Svg_SVG_reflection", "Svg_SVG__Body", "Drawing_DrawingOptions_reflection", "Drawing_DrawingOptions_New"]

