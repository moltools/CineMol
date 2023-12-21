from __future__ import annotations
from typing import Any
from .cine_mol_types import (Color, Chem_AtomType)
from .fable_modules.fable_library.double import divide
from .fable_modules.fable_library.reflection import (TypeInfo, union_type)
from .fable_modules.fable_library.types import (Array, Union)

def _expr1() -> TypeInfo:
    return union_type("CineMol.Style.AtomColorStyle", [], AtomColorStyle, lambda: [[]])


class AtomColorStyle(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["CPK"]


AtomColorStyle_reflection = _expr1

def AtomColorStyle__Color_Z3BD3A37D(this: AtomColorStyle, atom_type: Chem_AtomType) -> Color:
    tupled_arg: tuple[int, int, int] = ((255, 255, 255)) if (atom_type.tag == 0) else (((48, 48, 48)) if (atom_type.tag == 5) else (((0, 0, 255)) if (atom_type.tag == 6) else (((255, 0, 0)) if (atom_type.tag == 7) else (((255, 165, 0)) if (atom_type.tag == 14) else (((255, 255, 0)) if (atom_type.tag == 15) else (((245, 245, 220)) if (atom_type.tag == 4) else (((139, 0, 0)) if (atom_type.tag == 34) else (((148, 0, 211)) if (atom_type.tag == 52) else (((128, 128, 128)) if (atom_type.tag == 21) else (((255, 140, 0)) if (atom_type.tag == 25) else (((0, 128, 0)) if (atom_type.tag == 8) else (((0, 128, 0)) if (atom_type.tag == 16) else (((0, 255, 255)) if (atom_type.tag == 1) else (((0, 255, 255)) if (atom_type.tag == 9) else (((0, 255, 255)) if (atom_type.tag == 17) else (((0, 255, 255)) if (atom_type.tag == 35) else (((0, 255, 255)) if (atom_type.tag == 53) else (((238, 130, 238)) if (atom_type.tag == 2) else (((238, 130, 238)) if (atom_type.tag == 10) else (((238, 130, 238)) if (atom_type.tag == 18) else (((238, 130, 238)) if (atom_type.tag == 36) else (((238, 130, 238)) if (atom_type.tag == 54) else (((238, 130, 238)) if (atom_type.tag == 72) else (((0, 100, 0)) if (atom_type.tag == 3) else (((0, 100, 0)) if (atom_type.tag == 11) else (((0, 100, 0)) if (atom_type.tag == 19) else (((0, 100, 0)) if (atom_type.tag == 37) else (((0, 100, 0)) if (atom_type.tag == 55) else (((0, 100, 0)) if (atom_type.tag == 73) else ((255, 192, 203)))))))))))))))))))))))))))))))
    return Color(0, tupled_arg[0], tupled_arg[1], tupled_arg[2])


def _expr2() -> TypeInfo:
    return union_type("CineMol.Style.AtomRadius", [], AtomRadius, lambda: [[]])


class AtomRadius(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["PubChem"]


AtomRadius_reflection = _expr2

def AtomRadius__Radius_Z3BD3A37D(this: AtomRadius, atom_type: Chem_AtomType) -> float:
    return divide(140.0 if (atom_type.tag == 1) else (182.0 if (atom_type.tag == 2) else (153.0 if (atom_type.tag == 3) else (192.0 if (atom_type.tag == 4) else (170.0 if (atom_type.tag == 5) else (155.0 if (atom_type.tag == 6) else (152.0 if (atom_type.tag == 7) else (135.0 if (atom_type.tag == 8) else (154.0 if (atom_type.tag == 9) else (227.0 if (atom_type.tag == 10) else (173.0 if (atom_type.tag == 11) else (184.0 if (atom_type.tag == 12) else (210.0 if (atom_type.tag == 13) else (180.0 if (atom_type.tag == 14) else (180.0 if (atom_type.tag == 15) else (175.0 if (atom_type.tag == 16) else (188.0 if (atom_type.tag == 17) else (275.0 if (atom_type.tag == 18) else (231.0 if (atom_type.tag == 19) else (211.0 if (atom_type.tag == 20) else (187.0 if (atom_type.tag == 21) else (179.0 if (atom_type.tag == 22) else (189.0 if (atom_type.tag == 23) else (197.0 if (atom_type.tag == 24) else (194.0 if (atom_type.tag == 25) else (192.0 if (atom_type.tag == 26) else (163.0 if (atom_type.tag == 27) else (140.0 if (atom_type.tag == 28) else (139.0 if (atom_type.tag == 29) else (187.0 if (atom_type.tag == 30) else (211.0 if (atom_type.tag == 31) else (185.0 if (atom_type.tag == 32) else (190.0 if (atom_type.tag == 33) else (183.0 if (atom_type.tag == 34) else (202.0 if (atom_type.tag == 35) else (303.0 if (atom_type.tag == 36) else (249.0 if (atom_type.tag == 37) else (219.0 if (atom_type.tag == 38) else (186.0 if (atom_type.tag == 39) else (207.0 if (atom_type.tag == 40) else (209.0 if (atom_type.tag == 41) else (209.0 if (atom_type.tag == 42) else (207.0 if (atom_type.tag == 43) else (195.0 if (atom_type.tag == 44) else (202.0 if (atom_type.tag == 45) else (172.0 if (atom_type.tag == 46) else (158.0 if (atom_type.tag == 47) else (193.0 if (atom_type.tag == 48) else (217.0 if (atom_type.tag == 49) else (206.0 if (atom_type.tag == 50) else (206.0 if (atom_type.tag == 51) else (198.0 if (atom_type.tag == 52) else (216.0 if (atom_type.tag == 53) else (343.0 if (atom_type.tag == 54) else (268.0 if (atom_type.tag == 55) else (221.0 if (atom_type.tag == 56) else (212.0 if (atom_type.tag == 57) else (217.0 if (atom_type.tag == 58) else (210.0 if (atom_type.tag == 59) else (217.0 if (atom_type.tag == 60) else (216.0 if (atom_type.tag == 61) else (202.0 if (atom_type.tag == 62) else (209.0 if (atom_type.tag == 63) else (166.0 if (atom_type.tag == 64) else (209.0 if (atom_type.tag == 65) else (196.0 if (atom_type.tag == 66) else (202.0 if (atom_type.tag == 67) else (207.0 if (atom_type.tag == 68) else (197.0 if (atom_type.tag == 69) else (202.0 if (atom_type.tag == 70) else (220.0 if (atom_type.tag == 71) else (348.0 if (atom_type.tag == 72) else (283.0 if (atom_type.tag == 73) else 120.0)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))), 100.0)


__all__ = ["AtomColorStyle_reflection", "AtomColorStyle__Color_Z3BD3A37D", "AtomRadius_reflection", "AtomRadius__Radius_Z3BD3A37D"]

