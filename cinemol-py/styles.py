from __future__ import annotations
from array import array
from dataclasses import dataclass
from typing import (Tuple, Any, List)
from .fable_modules.fable_library.double import divide
from .fable_modules.fable_library.reflection import (TypeInfo, int32_type, tuple_type, float64_type, record_type, union_type)
from .fable_modules.fable_library.types import (Record, Array, Union)

def _expr1() -> TypeInfo:
    return record_type("Client.CineMol.Styles.Color", [], Color, lambda: [("RGB", tuple_type(int32_type, int32_type, int32_type)), ("Alpha", float64_type)])


@dataclass(eq = False, repr = False)
class Color(Record):
    RGB: Tuple[int, int, int]
    Alpha: float

Color_reflection = _expr1

def Color__Diffuse_5E38073B(x: Color, alpha: float) -> Tuple[float, float, float]:
    pattern_input: Tuple[int, int, int] = x.RGB
    return (pattern_input[0] * alpha, pattern_input[1] * alpha, pattern_input[2] * alpha)


atom_color_gradient: Array[float] = array("d", [1.0, 0.89, 0.66, 0.34, 0.0])

def _expr3() -> TypeInfo:
    return union_type("Client.CineMol.Styles.AtomType", [], AtomType, lambda: [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []])


class AtomType(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Unknown"]


AtomType_reflection = _expr3

def _expr4() -> TypeInfo:
    return union_type("Client.CineMol.Styles.AtomColorStyle", [], AtomColorStyle, lambda: [[]])


class AtomColorStyle(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["CPK"]


AtomColorStyle_reflection = _expr4

def get_atom_color(style: AtomColorStyle, atom_type: AtomType) -> Color:
    return Color(((255, 255, 255)) if (atom_type.tag == 0) else (((255, 255, 192)) if (atom_type.tag == 1) else (((204, 128, 255)) if (atom_type.tag == 2) else (((194, 255, 0)) if (atom_type.tag == 3) else (((255, 181, 181)) if (atom_type.tag == 4) else (((144, 144, 144)) if (atom_type.tag == 5) else (((48, 80, 248)) if (atom_type.tag == 6) else (((255, 13, 13)) if (atom_type.tag == 7) else (((144, 224, 80)) if (atom_type.tag == 8) else (((179, 227, 245)) if (atom_type.tag == 9) else (((171, 92, 242)) if (atom_type.tag == 10) else (((138, 255, 0)) if (atom_type.tag == 11) else (((191, 166, 166)) if (atom_type.tag == 12) else (((240, 200, 160)) if (atom_type.tag == 13) else (((255, 128, 0)) if (atom_type.tag == 14) else (((255, 255, 48)) if (atom_type.tag == 15) else (((31, 240, 31)) if (atom_type.tag == 16) else (((128, 209, 227)) if (atom_type.tag == 17) else (((143, 64, 212)) if (atom_type.tag == 18) else (((61, 255, 0)) if (atom_type.tag == 19) else (((230, 230, 230)) if (atom_type.tag == 20) else (((191, 194, 199)) if (atom_type.tag == 21) else (((166, 166, 171)) if (atom_type.tag == 22) else (((138, 153, 199)) if (atom_type.tag == 23) else (((156, 122, 199)) if (atom_type.tag == 24) else (((224, 102, 51)) if (atom_type.tag == 25) else (((240, 144, 160)) if (atom_type.tag == 26) else (((80, 208, 80)) if (atom_type.tag == 27) else (((200, 128, 51)) if (atom_type.tag == 28) else (((125, 128, 176)) if (atom_type.tag == 29) else (((194, 143, 143)) if (atom_type.tag == 30) else (((102, 143, 143)) if (atom_type.tag == 31) else (((189, 128, 227)) if (atom_type.tag == 32) else (((255, 161, 0)) if (atom_type.tag == 33) else (((166, 41, 41)) if (atom_type.tag == 34) else (((92, 184, 209)) if (atom_type.tag == 35) else (((112, 46, 176)) if (atom_type.tag == 36) else (((0, 255, 0)) if (atom_type.tag == 37) else (((148, 255, 255)) if (atom_type.tag == 38) else (((148, 224, 224)) if (atom_type.tag == 39) else (((115, 194, 201)) if (atom_type.tag == 40) else (((84, 181, 181)) if (atom_type.tag == 41) else (((59, 158, 158)) if (atom_type.tag == 42) else (((36, 143, 143)) if (atom_type.tag == 43) else (((10, 125, 140)) if (atom_type.tag == 44) else (((0, 105, 133)) if (atom_type.tag == 45) else (((192, 192, 192)) if (atom_type.tag == 46) else (((255, 217, 143)) if (atom_type.tag == 47) else (((166, 117, 115)) if (atom_type.tag == 48) else (((102, 128, 128)) if (atom_type.tag == 49) else (((158, 99, 181)) if (atom_type.tag == 50) else (((212, 122, 0)) if (atom_type.tag == 51) else (((148, 0, 148)) if (atom_type.tag == 52) else (((66, 158, 176)) if (atom_type.tag == 53) else (((87, 23, 143)) if (atom_type.tag == 54) else (((0, 201, 0)) if (atom_type.tag == 55) else (((112, 212, 255)) if (atom_type.tag == 56) else (((255, 255, 199)) if (atom_type.tag == 57) else (((217, 255, 199)) if (atom_type.tag == 58) else (((199, 255, 199)) if (atom_type.tag == 59) else (((163, 255, 199)) if (atom_type.tag == 60) else (((143, 255, 199)) if (atom_type.tag == 61) else (((97, 255, 199)) if (atom_type.tag == 62) else (((69, 255, 199)) if (atom_type.tag == 63) else (((48, 255, 199)) if (atom_type.tag == 64) else (((31, 255, 199)) if (atom_type.tag == 65) else (((0, 255, 156)) if (atom_type.tag == 66) else (((0, 230, 117)) if (atom_type.tag == 67) else (((0, 212, 82)) if (atom_type.tag == 68) else (((0, 191, 56)) if (atom_type.tag == 69) else (((0, 171, 36)) if (atom_type.tag == 70) else (((77, 194, 255)) if (atom_type.tag == 71) else (((77, 166, 255)) if (atom_type.tag == 72) else (((33, 148, 214)) if (atom_type.tag == 73) else (((38, 125, 171)) if (atom_type.tag == 74) else (((38, 102, 150)) if (atom_type.tag == 75) else (((23, 84, 135)) if (atom_type.tag == 76) else (((208, 208, 224)) if (atom_type.tag == 77) else (((255, 209, 35)) if (atom_type.tag == 78) else (((184, 184, 208)) if (atom_type.tag == 79) else (((166, 84, 77)) if (atom_type.tag == 80) else (((87, 89, 97)) if (atom_type.tag == 81) else (((158, 79, 181)) if (atom_type.tag == 82) else (((171, 92, 0)) if (atom_type.tag == 83) else (((117, 79, 69)) if (atom_type.tag == 84) else (((66, 130, 150)) if (atom_type.tag == 85) else (((66, 0, 102)) if (atom_type.tag == 86) else (((0, 125, 0)) if (atom_type.tag == 87) else (((112, 171, 250)) if (atom_type.tag == 88) else (((0, 186, 255)) if (atom_type.tag == 89) else (((0, 161, 255)) if (atom_type.tag == 90) else (((0, 143, 255)) if (atom_type.tag == 91) else (((0, 128, 255)) if (atom_type.tag == 92) else (((0, 107, 255)) if (atom_type.tag == 93) else (((84, 92, 242)) if (atom_type.tag == 94) else (((120, 92, 227)) if (atom_type.tag == 95) else (((138, 79, 227)) if (atom_type.tag == 96) else (((161, 54, 212)) if (atom_type.tag == 97) else (((179, 31, 212)) if (atom_type.tag == 98) else (((179, 31, 186)) if (atom_type.tag == 99) else (((179, 13, 166)) if (atom_type.tag == 100) else (((189, 13, 135)) if (atom_type.tag == 101) else (((199, 0, 102)) if (atom_type.tag == 102) else (((204, 0, 89)) if (atom_type.tag == 103) else (((209, 0, 79)) if (atom_type.tag == 104) else (((217, 0, 69)) if (atom_type.tag == 105) else (((224, 0, 56)) if (atom_type.tag == 106) else (((230, 0, 46)) if (atom_type.tag == 107) else (((235, 0, 38)) if (atom_type.tag == 108) else ((0, 0, 0)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))), 1.0)


def _expr5() -> TypeInfo:
    return union_type("Client.CineMol.Styles.AtomGeomStyle", [], AtomGeomStyle, lambda: [[]])


class AtomGeomStyle(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["Default"]


AtomGeomStyle_reflection = _expr5

def get_atom_radius(style: AtomGeomStyle, atom_type: AtomType) -> float:
    return 37 if (atom_type.tag == 0) else (152 if (atom_type.tag == 2) else (112 if (atom_type.tag == 3) else (88 if (atom_type.tag == 4) else (77 if (atom_type.tag == 5) else (70 if (atom_type.tag == 6) else (66 if (atom_type.tag == 7) else (64 if (atom_type.tag == 8) else (186 if (atom_type.tag == 10) else (160 if (atom_type.tag == 11) else (143 if (atom_type.tag == 12) else (117 if (atom_type.tag == 13) else (110 if (atom_type.tag == 14) else (104 if (atom_type.tag == 15) else (99 if (atom_type.tag == 16) else (231 if (atom_type.tag == 18) else (197 if (atom_type.tag == 19) else (160 if (atom_type.tag == 20) else (146 if (atom_type.tag == 21) else (131 if (atom_type.tag == 22) else (125 if (atom_type.tag == 23) else (129 if (atom_type.tag == 24) else (126 if (atom_type.tag == 25) else (126 if (atom_type.tag == 26) else (124 if (atom_type.tag == 27) else (128 if (atom_type.tag == 28) else (133 if (atom_type.tag == 29) else (122 if (atom_type.tag == 30) else (122 if (atom_type.tag == 31) else (121 if (atom_type.tag == 32) else (117 if (atom_type.tag == 33) else (114 if (atom_type.tag == 34) else (241 if (atom_type.tag == 36) else (215 if (atom_type.tag == 37) else (180 if (atom_type.tag == 38) else (157 if (atom_type.tag == 39) else (143 if (atom_type.tag == 40) else (136 if (atom_type.tag == 41) else (130 if (atom_type.tag == 42) else (133 if (atom_type.tag == 43) else (134 if (atom_type.tag == 44) else (138 if (atom_type.tag == 45) else (144 if (atom_type.tag == 46) else (149 if (atom_type.tag == 47) else (162 if (atom_type.tag == 48) else (140 if (atom_type.tag == 49) else (141 if (atom_type.tag == 50) else (137 if (atom_type.tag == 51) else (133 if (atom_type.tag == 52) else (262 if (atom_type.tag == 54) else (217 if (atom_type.tag == 55) else (157 if (atom_type.tag == 71) else (143 if (atom_type.tag == 72) else (137 if (atom_type.tag == 73) else (137 if (atom_type.tag == 74) else (134 if (atom_type.tag == 75) else (135 if (atom_type.tag == 76) else (138 if (atom_type.tag == 77) else (144 if (atom_type.tag == 78) else (150 if (atom_type.tag == 79) else (171 if (atom_type.tag == 80) else (175 if (atom_type.tag == 81) else (146 if (atom_type.tag == 82) else (140 if (atom_type.tag == 83) else (140 if (atom_type.tag == 84) else 77))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


def normalize_radius(style: AtomGeomStyle, radius: float) -> float:
    return divide(radius, get_atom_radius(AtomGeomStyle(0), AtomType(5)))


__all__ = ["Color_reflection", "Color__Diffuse_5E38073B", "atom_color_gradient", "AtomType_reflection", "AtomColorStyle_reflection", "get_atom_color", "AtomGeomStyle_reflection", "get_atom_radius", "normalize_radius"]

