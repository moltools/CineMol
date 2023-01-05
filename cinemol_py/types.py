from __future__ import annotations
from dataclasses import dataclass
from math import (cos, sin, pow, inf, acos)
from typing import (Any, List, Tuple, Callable, Optional)
from .fable_modules.fable_library.double import (sqrt, divide)
from .fable_modules.fable_library.list import FSharpList
from .fable_modules.fable_library.reflection import (TypeInfo, float64_type, record_type, union_type, tuple_type, int32_type, option_type, list_type)
from .fable_modules.fable_library.types import (Record, Array, Union)
from .styles import (Color, AtomType, AtomType_reflection, Color_reflection)

def _expr6() -> TypeInfo:
    return record_type("Client.CineMol.Types.Zoom", [], Zoom, lambda: [("Ratio", float64_type)])


@dataclass(eq = False, repr = False)
class Zoom(Record):
    Ratio: float

Zoom_reflection = _expr6

def Zoom_get_init(__unit: None=None) -> Zoom:
    return Zoom(1.0)


def _expr7() -> TypeInfo:
    return record_type("Client.CineMol.Types.Rotation", [], Rotation, lambda: [("AxisX", float64_type), ("AxisY", float64_type), ("AxisZ", float64_type)])


@dataclass(eq = False, repr = False)
class Rotation(Record):
    AxisX: float
    AxisY: float
    AxisZ: float

Rotation_reflection = _expr7

def Rotation_get_init(__unit: None=None) -> Rotation:
    return Rotation(0.0, 0.0, 0.0)


def _expr8() -> TypeInfo:
    return union_type("Client.CineMol.Types.Axis", [], Axis, lambda: [[], [], []])


class Axis(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["X", "Y", "Z"]


Axis_reflection = _expr8

def _expr9() -> TypeInfo:
    return record_type("Client.CineMol.Types.Point2D", [], Point2D, lambda: [("X", float64_type), ("Y", float64_type)])


@dataclass(eq = False, repr = False)
class Point2D(Record):
    X: float
    Y: float

Point2D_reflection = _expr9

def _expr10() -> TypeInfo:
    return record_type("Client.CineMol.Types.Vector2D", [], Vector2D, lambda: [("X", float64_type), ("Y", float64_type)])


@dataclass(eq = False, repr = False)
class Vector2D(Record):
    X: float
    Y: float

Vector2D_reflection = _expr10

def _expr11() -> TypeInfo:
    return record_type("Client.CineMol.Types.Point3D", [], Point3D, lambda: [("X", float64_type), ("Y", float64_type), ("Z", float64_type)])


@dataclass(eq = False, repr = False)
class Point3D(Record):
    X: float
    Y: float
    Z: float

Point3D_reflection = _expr11

def _expr12() -> TypeInfo:
    return record_type("Client.CineMol.Types.Vector3D", [], Vector3D, lambda: [("X", float64_type), ("Y", float64_type), ("Z", float64_type)])


@dataclass(eq = False, repr = False)
class Vector3D(Record):
    X: float
    Y: float
    Z: float

Vector3D_reflection = _expr12

def Axis__get_RotationMatrix(x: Axis) -> Callable[[Tuple[Point3D, float]], Point3D]:
    if x.tag == 1:
        def _arrow13(tupled_arg_1: Tuple[Point3D, float], x: Axis=x) -> Point3D:
            p_1: Point3D = tupled_arg_1[0]
            rad_1: float = tupled_arg_1[1]
            return Point3D((p_1.X * cos(rad_1)) + (p_1.Z * sin(rad_1)), p_1.Y, (p_1.Z * cos(rad_1)) - (p_1.X * sin(rad_1)))

        return _arrow13

    elif x.tag == 2:
        def _arrow14(tupled_arg_2: Tuple[Point3D, float], x: Axis=x) -> Point3D:
            p_2: Point3D = tupled_arg_2[0]
            rad_2: float = tupled_arg_2[1]
            return Point3D((p_2.X * cos(rad_2)) - (p_2.Y * sin(rad_2)), (p_2.X * sin(rad_2)) + (p_2.Y * cos(rad_2)), p_2.Z)

        return _arrow14

    else: 
        def _arrow15(tupled_arg: Tuple[Point3D, float], x: Axis=x) -> Point3D:
            p: Point3D = tupled_arg[0]
            rad: float = tupled_arg[1]
            return Point3D(p.X, (p.Y * cos(rad)) - (p.Z * sin(rad)), (p.Y * sin(rad)) + (p.Z * cos(rad)))

        return _arrow15



def Point2D_op_Subtraction_25FD1980(p1: Point2D, p2: Point2D) -> Point2D:
    return Point2D(p1.X - p2.X, p1.Y - p2.Y)


def Point2D_Pow(p: Point2D, d: float) -> Point2D:
    return Point2D(pow(p.X, d), pow(p.Y, d))


def Point2D_Sum_591E284C(p: Point2D) -> float:
    return p.X + p.Y


def Point2D__Distance_591E284C(p1: Point2D, p2: Point2D) -> float:
    return sqrt(Point2D_Sum_591E284C(Point2D_Pow(Point2D_op_Subtraction_25FD1980(p1, p2), 2.0)))


def Point2D__Midpoint_591E284C(p1: Point2D, p2: Point2D) -> Point2D:
    return Point2D(divide(p1.X + p2.X, 2.0), divide(p1.Y + p2.Y, 2.0))


def Point2D__FindVector_591E284C(p1: Point2D, p2: Point2D) -> Vector2D:
    return Vector2D(p2.X - p1.X, p2.Y - p1.Y)


def Vector2D__get_SumOfSquares(u: Vector2D) -> float:
    return pow(u.X, 2.0) + pow(u.Y, 2.0)


def Vector2D__get_Magnitude(u: Vector2D) -> float:
    return sqrt(Vector2D__get_SumOfSquares(u))


def Vector2D_op_Multiply_69088842(k: float, v: Vector2D) -> Vector2D:
    return Vector2D(k * v.X, k * v.Y)


def Vector2D_op_Addition_Z61F2D8E0(v1: Vector2D, v2: Vector2D) -> Vector2D:
    return Vector2D(v1.X + v2.X, v1.Y + v2.Y)


def Vector2D__get_Norm(u: Vector2D) -> Vector2D:
    mag: float = Vector2D__get_Magnitude(u)
    return Vector2D_op_Multiply_69088842(inf if (mag == 0.0) else divide(1.0, mag), u)


def Vector2D__Dot_4C3066D9(u: Vector2D, v: Vector2D) -> float:
    return (u.X * v.X) + (u.Y * v.Y)


def Point3D_op_Subtraction_25FD1E60(p1: Point3D, p2: Point3D) -> Point3D:
    return Point3D(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z)


def Point3D_Pow(p: Point3D, d: float) -> Point3D:
    return Point3D(pow(p.X, d), pow(p.Y, d), pow(p.Z, d))


def Point3D_Sum_591E286D(p: Point3D) -> float:
    return (p.X + p.Y) + p.Z


def Point3D__Distance_591E286D(p1: Point3D, p2: Point3D) -> float:
    return sqrt(Point3D_Sum_591E286D(Point3D_Pow(Point3D_op_Subtraction_25FD1E60(p1, p2), 2.0)))


def Point3D__Centroid_591E286D(p1: Point3D, p2: Point3D) -> Point3D:
    return Point3D(divide(p1.X + p2.X, 2.0), divide(p1.Y + p2.X, 2.0), divide(p1.Z + p2.Z, 2.0))


def Point3D__Rotate(p: Point3D, axis: Axis, rad: float) -> Point3D:
    return Axis__get_RotationMatrix(axis)((p, rad))


def Point3D__FindVector_591E286D(p1: Point3D, p2: Point3D) -> Vector3D:
    return Vector3D(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z)


def Vector3D__get_SumOfSquares(u: Vector3D) -> float:
    return (pow(u.X, 2.0) + pow(u.Y, 2.0)) + pow(u.Z, 2.0)


def Vector3D__get_Magnitude(u: Vector3D) -> float:
    return sqrt(Vector3D__get_SumOfSquares(u))


def Vector3D_op_Multiply_690888A3(k: float, v: Vector3D) -> Vector3D:
    return Vector3D(k * v.X, k * v.Y, k * v.Z)


def Vector3D_op_Addition_Z61F2B500(v1: Vector3D, v2: Vector3D) -> Vector3D:
    return Vector3D(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z)


def Vector3D__get_Norm(u: Vector3D) -> Vector3D:
    mag: float = Vector3D__get_Magnitude(u)
    return Vector3D_op_Multiply_690888A3(inf if (mag == 0.0) else divide(1.0, mag), u)


def Vector3D__Dot_4C306638(u: Vector3D, v: Vector3D) -> float:
    return ((u.X * v.X) + (u.Y * v.Y)) + (u.Z * v.Z)


def Vector3D__Cross_4C306638(u: Vector3D, v: Vector3D) -> Vector3D:
    return Vector3D((u.Y * v.Z) - (u.Z * v.Y), (u.Z * v.X) - (u.X * v.Z), (u.X * v.Y) - (u.Y * v.X))


def Vector3D__ProjectVector_4C306638(x: Vector3D, v: Vector3D) -> float:
    return divide(Vector3D__Dot_4C306638(v, x), Vector3D__get_Magnitude(v))


def _expr16() -> TypeInfo:
    return union_type("Client.CineMol.Types.SphereSphereIntersection", [], SphereSphereIntersection, lambda: [[], [], [("Item", Point3D_reflection())], [("Item1", Point3D_reflection()), ("Item2", float64_type), ("Item3", Vector3D_reflection())]])


class SphereSphereIntersection(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["Eclipsed", "NoIntersection", "IntersectionPoint", "IntersectionCircle"]


SphereSphereIntersection_reflection = _expr16

def _expr17() -> TypeInfo:
    return record_type("Client.CineMol.Types.ClipPath", [], ClipPath, lambda: [("Line", tuple_type(Point2D_reflection(), Point2D_reflection())), ("SelectForSide", SelectForSide_reflection())])


@dataclass(eq = False, repr = False)
class ClipPath(Record):
    Line: Tuple[Point2D, Point2D]
    SelectForSide: SelectForSide

ClipPath_reflection = _expr17

def _expr18() -> TypeInfo:
    return union_type("Client.CineMol.Types.SelectForSide", [], SelectForSide, lambda: [[("Item", Point2D_reflection())], [("Item", Point2D_reflection())], []])


class SelectForSide(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["IncludeSide", "ExcludeSide", "IncludeBothSides"]


SelectForSide_reflection = _expr18

def _expr19() -> TypeInfo:
    return record_type("Client.CineMol.Types.AtomInfo", [], AtomInfo, lambda: [("Index", int32_type), ("AtomType", AtomType_reflection()), ("Center", Point3D_reflection()), ("Radius", float64_type), ("Color", option_type(Color_reflection()))])


@dataclass(eq = False, repr = False)
class AtomInfo(Record):
    Index: int
    AtomType: AtomType
    Center: Point3D
    Radius: float
    Color: Optional[Color]

AtomInfo_reflection = _expr19

def AtomInfo__Rotate(x: AtomInfo, axis: Axis, rad: float) -> AtomInfo:
    return AtomInfo(x.Index, x.AtomType, Point3D__Rotate(x.Center, axis, rad), x.Radius, x.Color)


def AtomInfo__Intersects_Z6467FD91(this: AtomInfo, other: AtomInfo) -> bool:
    if Point3D__Distance_591E286D(this.Center, other.Center) <= (this.Radius + other.Radius):
        return True

    else: 
        return False



def AtomInfo__Intersection_Z6467FD91(this: AtomInfo, other: AtomInfo) -> SphereSphereIntersection:
    dist: float = Point3D__Distance_591E286D(this.Center, other.Center)
    def _arrow20(__unit: None=None, this: AtomInfo=this, other: AtomInfo=other) -> bool:
        d: float = dist
        return True if (d >= (this.Radius + other.Radius)) else ((this.Radius == other.Radius) if (d == 0.0) else False)

    if _arrow20():
        return SphereSphereIntersection(1)

    elif (dist + this.Radius) < other.Radius:
        return SphereSphereIntersection(0)

    else: 
        A: float = 2.0 * (other.Center.X - this.Center.X)
        B: float = 2.0 * (other.Center.Y - this.Center.Y)
        C: float = 2.0 * (other.Center.Z - this.Center.Z)
        t: float = divide((((this.Center.X * A) + (this.Center.Y * B)) + (this.Center.Z * C)) + (((((((pow(this.Center.X, 2.0) - pow(other.Center.X, 2.0)) + pow(this.Center.Y, 2.0)) - pow(other.Center.Y, 2.0)) + pow(this.Center.Z, 2.0)) - pow(other.Center.Z, 2.0)) - pow(this.Radius, 2.0)) + pow(other.Radius, 2.0)), ((A * (this.Center.X - other.Center.X)) + (B * (this.Center.Y - other.Center.Y))) + (C * (this.Center.Z - other.Center.Z)))
        intersection_center: Point3D = Point3D(this.Center.X + (t * (other.Center.X - this.Center.X)), this.Center.Y + (t * (other.Center.Y - this.Center.Y)), this.Center.Z + (t * (other.Center.Z - this.Center.Z)))
        x_1: float = divide((pow(this.Radius, 2.0) + pow(dist, 2.0)) - pow(other.Radius, 2.0), (2.0 * this.Radius) * dist)
        if x_1 < 1.0:
            R: float = this.Radius * sin(acos(x_1))
            if R == 0.0:
                return SphereSphereIntersection(2, intersection_center)

            else: 
                return SphereSphereIntersection(3, intersection_center, R, Point3D__FindVector_591E286D(this.Center, other.Center))


        else: 
            return SphereSphereIntersection(1)




def _expr21() -> TypeInfo:
    return record_type("Client.CineMol.Types.ProjectedAtomInfo", [], ProjectedAtomInfo, lambda: [("Index", int32_type), ("AtomType", AtomType_reflection()), ("Center", Point2D_reflection()), ("Radius", float64_type), ("ClipPaths", list_type(ClipPath_reflection())), ("Color", option_type(Color_reflection()))])


@dataclass(eq = False, repr = False)
class ProjectedAtomInfo(Record):
    Index: int
    AtomType: AtomType
    Center: Point2D
    Radius: float
    ClipPaths: FSharpList[ClipPath]
    Color: Optional[Color]

ProjectedAtomInfo_reflection = _expr21

def _expr22() -> TypeInfo:
    return record_type("Client.CineMol.Types.BondInfo", [], BondInfo, lambda: [("Index", int32_type), ("Start", int32_type), ("End", int32_type), ("BondType", BondType_reflection()), ("Scaling", float64_type), ("StartColor", option_type(Color_reflection())), ("EndColor", option_type(Color_reflection()))])


@dataclass(eq = False, repr = False)
class BondInfo(Record):
    Index: int
    Start: int
    End: int
    BondType: BondType
    Scaling: float
    StartColor: Optional[Color]
    EndColor: Optional[Color]

BondInfo_reflection = _expr22

def _expr23() -> TypeInfo:
    return union_type("Client.CineMol.Types.BondType", [], BondType, lambda: [[], [], [], [], []])


class BondType(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["Single", "Double", "Triple", "Aromatic", "Unknown"]


BondType_reflection = _expr23

def create_atom(index: int, atom_type: AtomType, c: Point3D, r: float) -> AtomInfo:
    return AtomInfo(index, atom_type, c, r, None)


def _expr24() -> TypeInfo:
    return record_type("Client.CineMol.Types.Molecule", [], Molecule, lambda: [("Atoms", list_type(AtomInfo_reflection())), ("Bonds", list_type(BondInfo_reflection()))])


@dataclass(eq = False, repr = False)
class Molecule(Record):
    Atoms: FSharpList[AtomInfo]
    Bonds: FSharpList[BondInfo]

Molecule_reflection = _expr24

def _expr25() -> TypeInfo:
    return record_type("Client.CineMol.Types.ProjectedMolecule", [], ProjectedMolecule, lambda: [("Atoms", list_type(ProjectedAtomInfo_reflection())), ("Bonds", list_type(BondInfo_reflection()))])


@dataclass(eq = False, repr = False)
class ProjectedMolecule(Record):
    Atoms: FSharpList[ProjectedAtomInfo]
    Bonds: FSharpList[BondInfo]

ProjectedMolecule_reflection = _expr25

def _expr26() -> TypeInfo:
    return union_type("Client.CineMol.Types.Depiction", [], Depiction, lambda: [[], [], []])


class Depiction(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["Filled", "BallAndStick", "Wire"]


Depiction_reflection = _expr26

origin: Point3D = Point3D(0.0, 0.0, 0.0)

def physical_projection(camera_perpendicular: Vector3D, camera_horizon: Vector3D, camera_forward: Vector3D, pov: Point3D, p: Point3D) -> Point3D:
    point_vector: Vector3D = Point3D__FindVector_591E286D(pov, p)
    return Point3D(Vector3D__ProjectVector_4C306638(point_vector, camera_perpendicular), Vector3D__ProjectVector_4C306638(point_vector, camera_horizon), Vector3D__ProjectVector_4C306638(point_vector, camera_forward))


def perspective_projection(focal_length: float, p: Point3D) -> Point3D:
    scale_factor: float = divide(focal_length, p.Z)
    return Point3D(p.X * scale_factor, p.Y * scale_factor, p.Z * scale_factor)


def project(camera_perpendicular: Vector3D, camera_horizon: Vector3D, camera_forward: Vector3D, pov: Point3D, focal_length: float, p: Point3D) -> Point3D:
    return perspective_projection(focal_length, physical_projection(camera_perpendicular, camera_horizon, camera_forward, pov, p))


__all__ = ["Zoom_reflection", "Zoom_get_init", "Rotation_reflection", "Rotation_get_init", "Axis_reflection", "Point2D_reflection", "Vector2D_reflection", "Point3D_reflection", "Vector3D_reflection", "Axis__get_RotationMatrix", "Point2D_op_Subtraction_25FD1980", "Point2D_Pow", "Point2D_Sum_591E284C", "Point2D__Distance_591E284C", "Point2D__Midpoint_591E284C", "Point2D__FindVector_591E284C", "Vector2D__get_SumOfSquares", "Vector2D__get_Magnitude", "Vector2D_op_Multiply_69088842", "Vector2D_op_Addition_Z61F2D8E0", "Vector2D__get_Norm", "Vector2D__Dot_4C3066D9", "Point3D_op_Subtraction_25FD1E60", "Point3D_Pow", "Point3D_Sum_591E286D", "Point3D__Distance_591E286D", "Point3D__Centroid_591E286D", "Point3D__Rotate", "Point3D__FindVector_591E286D", "Vector3D__get_SumOfSquares", "Vector3D__get_Magnitude", "Vector3D_op_Multiply_690888A3", "Vector3D_op_Addition_Z61F2B500", "Vector3D__get_Norm", "Vector3D__Dot_4C306638", "Vector3D__Cross_4C306638", "Vector3D__ProjectVector_4C306638", "SphereSphereIntersection_reflection", "ClipPath_reflection", "SelectForSide_reflection", "AtomInfo_reflection", "AtomInfo__Rotate", "AtomInfo__Intersects_Z6467FD91", "AtomInfo__Intersection_Z6467FD91", "ProjectedAtomInfo_reflection", "BondInfo_reflection", "BondType_reflection", "create_atom", "Molecule_reflection", "ProjectedMolecule_reflection", "Depiction_reflection", "origin", "physical_projection", "perspective_projection", "project"]

