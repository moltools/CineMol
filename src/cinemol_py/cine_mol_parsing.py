from __future__ import annotations
from typing import Any
from .cine_mol_style import (AtomColorStyle__Color_Z3BD3A37D, AtomColorStyle, AtomRadius__Radius_Z3BD3A37D, AtomRadius)
from .cine_mol_types import (Chem_Atom, Chem_Bond, Chem_Molecule__AdjustForCentroid, Chem_Molecule, Chem_AtomType_FromString_Z721C83C5, Chem_AtomType, Geometry_Point3D, Chem_BondType_FromString_Z721C83C5, Chem_BondType)
from .fable_modules.fable_library.double import parse
from .fable_modules.fable_library.int32 import parse as parse_1
from .fable_modules.fable_library.list import (append, singleton, FSharpList, empty, map, length, tail, is_empty, head)
from .fable_modules.fable_library.range import range_big_int
from .fable_modules.fable_library.reflection import (TypeInfo, class_type, union_type)
from .fable_modules.fable_library.reg_exp import (match, create, groups)
from .fable_modules.fable_library.seq import (to_list, delay, collect, singleton as singleton_1, append as append_1, empty as empty_1, map as map_1)
from .fable_modules.fable_library.string_ import join
from .fable_modules.fable_library.types import (FSharpException, Array, Union)
from .fable_modules.fable_library.util import IEnumerable_1

def _expr23() -> TypeInfo:
    return class_type("CineMol.Parsing.ParserError", None, ParserError, class_type("System.Exception"))


class ParserError(FSharpException):
    def __init__(self, Data0: str) -> None:
        super().__init__()
        self.Data0 = Data0


ParserError_reflection = _expr23

def _expr24() -> TypeInfo:
    return union_type("CineMol.Parsing.FileType", [], FileType, lambda: [[]])


class FileType(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["Sdf"]


FileType_reflection = _expr24

def _expr25() -> TypeInfo:
    return union_type("CineMol.Parsing.FileParser", [], FileParser, lambda: [[("Item", FileType_reflection())]])


class FileParser(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> list[str]:
        return ["FileParser"]


FileParser_reflection = _expr25

def FileParser__Parse_Z721C83C5(this: FileParser, file_content: str) -> FSharpList[Chem_Molecule] | None:
    d_cap: str = ("(" + "[-+]?[0-9]*\\.?[0-9]+") + ")"
    def _arrow34(__unit: None=None, this: Any=this, file_content: Any=file_content) -> FSharpList[str]:
        def _arrow31(__unit: None=None) -> IEnumerable_1[str]:
            def _arrow30(match_value_1: int) -> IEnumerable_1[str]:
                return singleton_1("\\s{1,}[-+]?[0-9]*\\.?[0-9]+")

            return collect(_arrow30, to_list(range_big_int(0, 1, 11)))

        list_6: FSharpList[str] = append(singleton("\\s{1,}(\\w+)"), append(to_list(delay(_arrow31)), singleton("$")))
        def _arrow33(__unit: None=None) -> IEnumerable_1[str]:
            def _arrow32(match_value_2: int) -> IEnumerable_1[str]:
                return singleton_1("\\s{1,}" + d_cap)

            return collect(_arrow32, to_list(range_big_int(0, 1, 2)))

        return append(to_list(delay(_arrow33)), list_6)

    atom_line: str = join("", append(singleton("^"), _arrow34()))
    def _arrow39(__unit: None=None, this: Any=this, file_content: Any=file_content) -> FSharpList[str]:
        def _arrow36(__unit: None=None) -> IEnumerable_1[str]:
            def _arrow35(match_value_3: int) -> IEnumerable_1[str]:
                return singleton_1("\\s{1,}[-+]?[0-9]*\\.?[0-9]+")

            return collect(_arrow35, to_list(range_big_int(0, 1, 3)))

        list_12: FSharpList[str] = append(to_list(delay(_arrow36)), singleton("$"))
        def _arrow38(__unit: None=None) -> IEnumerable_1[str]:
            def _arrow37(match_value_4: int) -> IEnumerable_1[str]:
                return singleton_1("\\s{1,}" + d_cap)

            return collect(_arrow37, to_list(range_big_int(0, 1, 2)))

        return append(to_list(delay(_arrow38)), list_12)

    bond_line: str = join("", append(singleton("^"), _arrow39()))
    count: int = 0
    atoms: FSharpList[Chem_Atom] = empty()
    bonds: FSharpList[Chem_Bond] = empty()
    def mapping(mol: Chem_Molecule, this: Any=this, file_content: Any=file_content) -> Chem_Molecule:
        return Chem_Molecule__AdjustForCentroid(mol)

    def _arrow47(__unit: None=None, this: Any=this, file_content: Any=file_content) -> IEnumerable_1[Chem_Molecule]:
        def _arrow46(line_2: str) -> IEnumerable_1[Chem_Molecule]:
            nonlocal atoms, bonds
            match_value_5: str = line_2
            (pattern_matching_result, atom_symbol, x, y, z) = (None, None, None, None, None)
            if (match_value_5.find("$$$$") >= 0) == True:
                pattern_matching_result = 0

            else: 
                active_pattern_result: FSharpList[str] | None
                m: Any = match(create(atom_line), match_value_5)
                def _arrow45(__unit: None=None) -> IEnumerable_1[str]:
                    def _arrow44(g: Any) -> str:
                        return g or ""

                    return map_1(_arrow44, groups(m))

                active_pattern_result = tail(to_list(delay(_arrow45))) if (m is not None) else None
                if active_pattern_result is not None:
                    if not is_empty(active_pattern_result):
                        if not is_empty(tail(active_pattern_result)):
                            if not is_empty(tail(tail(active_pattern_result))):
                                if not is_empty(tail(tail(tail(active_pattern_result)))):
                                    if is_empty(tail(tail(tail(tail(active_pattern_result))))):
                                        pattern_matching_result = 1
                                        atom_symbol = head(tail(tail(tail(active_pattern_result))))
                                        x = head(active_pattern_result)
                                        y = head(tail(active_pattern_result))
                                        z = head(tail(tail(active_pattern_result)))

                                    else: 
                                        pattern_matching_result = 2


                                else: 
                                    pattern_matching_result = 2


                            else: 
                                pattern_matching_result = 2


                        else: 
                            pattern_matching_result = 2


                    else: 
                        pattern_matching_result = 2


                else: 
                    pattern_matching_result = 2


            if pattern_matching_result == 0:
                def _arrow41(__unit: None=None) -> IEnumerable_1[Chem_Molecule]:
                    nonlocal count, atoms, bonds
                    count = (count + 1) or 0
                    atoms = empty()
                    bonds = empty()
                    return empty_1()

                return append_1(singleton_1(Chem_Molecule(atoms, bonds)), delay(_arrow41))

            elif pattern_matching_result == 1:
                matchValue: float | None
                try: 
                    matchValue = parse(x)

                except Exception as match_value_6:
                    matchValue = None

                matchValue_1: float | None
                try: 
                    matchValue_1 = parse(y)

                except Exception as match_value_7:
                    matchValue_1 = None

                matchValue_2: float | None
                try: 
                    matchValue_2 = parse(z)

                except Exception as match_value_8:
                    matchValue_2 = None

                matchValue_3: Chem_AtomType | None = Chem_AtomType_FromString_Z721C83C5(atom_symbol)
                (pattern_matching_result_1, atom_type, x_1, y_1, z_1) = (None, None, None, None, None)
                if matchValue is not None:
                    if matchValue_1 is not None:
                        if matchValue_2 is not None:
                            if matchValue_3 is not None:
                                pattern_matching_result_1 = 0
                                atom_type = matchValue_3
                                x_1 = matchValue
                                y_1 = matchValue_1
                                z_1 = matchValue_2

                            else: 
                                pattern_matching_result_1 = 1


                        else: 
                            pattern_matching_result_1 = 1


                    else: 
                        pattern_matching_result_1 = 1


                else: 
                    pattern_matching_result_1 = 1

                if pattern_matching_result_1 == 0:
                    atom: Chem_Atom = Chem_Atom((length(atoms) + (length(bonds) * 2)) + 1, atom_type, AtomColorStyle__Color_Z3BD3A37D(AtomColorStyle(0), atom_type), 1.0, Geometry_Point3D(x_1, y_1, z_1), AtomRadius__Radius_Z3BD3A37D(AtomRadius(0), atom_type))
                    atoms = append(atoms, singleton(atom))
                    return empty_1()

                elif pattern_matching_result_1 == 1:
                    raise ParserError(("Unable to parse atom line: \'" + line_2) + "\'")
                    return empty_1()


            elif pattern_matching_result == 2:
                (pattern_matching_result_2, bond_type, e_idx, s_idx) = (None, None, None, None)
                active_pattern_result_1: FSharpList[str] | None
                m_1: Any = match(create(bond_line), match_value_5)
                def _arrow43(__unit: None=None) -> IEnumerable_1[str]:
                    def _arrow42(g_1: Any) -> str:
                        return g_1 or ""

                    return map_1(_arrow42, groups(m_1))

                active_pattern_result_1 = tail(to_list(delay(_arrow43))) if (m_1 is not None) else None
                if active_pattern_result_1 is not None:
                    if not is_empty(active_pattern_result_1):
                        if not is_empty(tail(active_pattern_result_1)):
                            if not is_empty(tail(tail(active_pattern_result_1))):
                                if is_empty(tail(tail(tail(active_pattern_result_1)))):
                                    pattern_matching_result_2 = 0
                                    bond_type = head(tail(tail(active_pattern_result_1)))
                                    e_idx = head(tail(active_pattern_result_1))
                                    s_idx = head(active_pattern_result_1)

                                else: 
                                    pattern_matching_result_2 = 1


                            else: 
                                pattern_matching_result_2 = 1


                        else: 
                            pattern_matching_result_2 = 1


                    else: 
                        pattern_matching_result_2 = 1


                else: 
                    pattern_matching_result_2 = 1

                if pattern_matching_result_2 == 0:
                    matchValue_4: int | None
                    try: 
                        matchValue_4 = parse_1(s_idx, 511, False, 32)

                    except Exception as match_value_10:
                        matchValue_4 = None

                    matchValue_5: int | None
                    try: 
                        matchValue_5 = parse_1(e_idx, 511, False, 32)

                    except Exception as match_value_11:
                        matchValue_5 = None

                    matchValue_6: Chem_BondType | None = Chem_BondType_FromString_Z721C83C5(bond_type)
                    (pattern_matching_result_3, bond_type_1, e_idx_1, s_idx_1) = (None, None, None, None)
                    if matchValue_4 is not None:
                        if matchValue_5 is not None:
                            if matchValue_6 is not None:
                                pattern_matching_result_3 = 0
                                bond_type_1 = matchValue_6
                                e_idx_1 = matchValue_5
                                s_idx_1 = matchValue_4

                            else: 
                                pattern_matching_result_3 = 1


                        else: 
                            pattern_matching_result_3 = 1


                    else: 
                        pattern_matching_result_3 = 1

                    if pattern_matching_result_3 == 0:
                        bond: Chem_Bond = Chem_Bond((length(atoms) + (length(bonds) * 2)) + 1, (length(atoms) + (length(bonds) * 2)) + 2, bond_type_1, s_idx_1, e_idx_1, None, None, 0.5)
                        bonds = append(bonds, singleton(bond))
                        return empty_1()

                    elif pattern_matching_result_3 == 1:
                        raise ParserError(("Unable to parse bond line: \'" + line_2) + "\'")
                        return empty_1()


                elif pattern_matching_result_2 == 1:
                    return empty_1()



        return collect(_arrow46, file_content.split("\n"))

    return map(mapping, to_list(delay(_arrow47)))


__all__ = ["ParserError_reflection", "FileType_reflection", "FileParser_reflection", "FileParser__Parse_Z721C83C5"]

