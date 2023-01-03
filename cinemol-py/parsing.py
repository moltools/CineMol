from __future__ import annotations
from typing import (Any, Optional)
from .fable_modules.fable_library.double import parse
from .fable_modules.fable_library.int32 import parse as parse_1
from .fable_modules.fable_library.list import (append, singleton, FSharpList, tail, empty, of_array, is_empty, head)
from .fable_modules.fable_library.range import range_big_int
from .fable_modules.fable_library.reflection import (TypeInfo, class_type)
from .fable_modules.fable_library.reg_exp import (match, groups)
from .fable_modules.fable_library.seq import (to_list, delay, collect, singleton as singleton_1, map, to_array, append as append_1, empty as empty_1)
from .fable_modules.fable_library.string import join
from .fable_modules.fable_library.types import (FSharpException, Array)
from .fable_modules.fable_library.util import (equals, IEnumerable_1)
from .styles import (AtomType, normalize_radius, AtomGeomStyle, get_atom_radius)
from .types import (BondType, AtomInfo, BondInfo, Molecule, Point3D, create_atom)

def _expr27() -> TypeInfo:
    return class_type("Client.CineMol.Parsing.ParserError", None, ParserError, class_type("System.Exception"))


class ParserError(FSharpException):
    def __init__(self, Data0: str) -> None:
        super().__init__()
        self.Data0 = Data0


ParserError_reflection = _expr27

def ParserError__Equals_229D3F39(this: Exception, obj: Exception) -> bool:
    if not equals(this, None):
        if not equals(obj, None):
            if isinstance(obj, ParserError):
                return this.Data0 == obj.Data0

            else: 
                return False


        else: 
            return False


    elif not equals(obj, None):
        return False

    else: 
        return True



s: str = "\\s{1,}"

d: str = "[-+]?[0-9]*\\.?[0-9]+"

d_cap: str = ("(" + d) + ")"

w_cap: str = "(\\w+)"

def _arrow32(__unit: None=None) -> FSharpList[str]:
    def _arrow29(__unit: None=None) -> IEnumerable_1[str]:
        def _arrow28(match_value: int) -> IEnumerable_1[str]:
            return singleton_1(s + d)

        return collect(_arrow28, to_list(range_big_int(0, 1, 11)))

    list_6: FSharpList[str] = append(singleton(s + w_cap), append(to_list(delay(_arrow29)), singleton("$")))
    def _arrow31(__unit: None=None) -> IEnumerable_1[str]:
        def _arrow30(match_value_1: int) -> IEnumerable_1[str]:
            return singleton_1(s + d_cap)

        return collect(_arrow30, to_list(range_big_int(0, 1, 2)))

    return append(to_list(delay(_arrow31)), list_6)


atom_line: str = join("", append(singleton("^"), _arrow32()))

def _007CAtomLine_007C__007C(line: str) -> Optional[FSharpList[str]]:
    m: Any = match(line, atom_line)
    if m is not None:
        def _arrow34(__unit: None=None, line: str=line) -> IEnumerable_1[str]:
            def _arrow33(g: Any) -> str:
                return g or ""

            return map(_arrow33, groups(m))

        return tail(to_list(delay(_arrow34)))

    else: 
        return None



def _arrow39(__unit: None=None) -> FSharpList[str]:
    def _arrow36(__unit: None=None) -> IEnumerable_1[str]:
        def _arrow35(match_value: int) -> IEnumerable_1[str]:
            return singleton_1(s + d)

        return collect(_arrow35, to_list(range_big_int(0, 1, 3)))

    list_4: FSharpList[str] = append(to_list(delay(_arrow36)), singleton("$"))
    def _arrow38(__unit: None=None) -> IEnumerable_1[str]:
        def _arrow37(match_value_1: int) -> IEnumerable_1[str]:
            return singleton_1(s + d_cap)

        return collect(_arrow37, to_list(range_big_int(0, 1, 2)))

    return append(to_list(delay(_arrow38)), list_4)


bond_line: str = join("", append(singleton("^"), _arrow39()))

def _007CBondLine_007C__007C(line: str) -> Optional[FSharpList[str]]:
    m: Any = match(line, bond_line)
    if m is not None:
        def _arrow41(__unit: None=None, line: str=line) -> IEnumerable_1[str]:
            def _arrow40(g: Any) -> str:
                return g or ""

            return map(_arrow40, groups(m))

        return tail(to_list(delay(_arrow41)))

    else: 
        return None



def try_cast_to_atom(atom: str) -> AtomType:
    if atom == "H":
        return AtomType(0)

    elif atom == "He":
        return AtomType(1)

    elif atom == "Li":
        return AtomType(2)

    elif atom == "Be":
        return AtomType(3)

    elif atom == "B":
        return AtomType(4)

    elif atom == "C":
        return AtomType(5)

    elif atom == "N":
        return AtomType(6)

    elif atom == "O":
        return AtomType(7)

    elif atom == "F":
        return AtomType(8)

    elif atom == "Ne":
        return AtomType(9)

    elif atom == "Na":
        return AtomType(10)

    elif atom == "Mg":
        return AtomType(11)

    elif atom == "Al":
        return AtomType(12)

    elif atom == "Si":
        return AtomType(13)

    elif atom == "P":
        return AtomType(14)

    elif atom == "S":
        return AtomType(15)

    elif atom == "Cl":
        return AtomType(16)

    elif atom == "Ar":
        return AtomType(17)

    elif atom == "K":
        return AtomType(18)

    elif atom == "Ca":
        return AtomType(19)

    elif atom == "Sc":
        return AtomType(20)

    elif atom == "Ti":
        return AtomType(21)

    elif atom == "V":
        return AtomType(22)

    elif atom == "Cr":
        return AtomType(23)

    elif atom == "Mn":
        return AtomType(24)

    elif atom == "Fe":
        return AtomType(25)

    elif atom == "Co":
        return AtomType(26)

    elif atom == "Ni":
        return AtomType(27)

    elif atom == "Cu":
        return AtomType(28)

    elif atom == "Zn":
        return AtomType(29)

    elif atom == "Ga":
        return AtomType(30)

    elif atom == "Ge":
        return AtomType(31)

    elif atom == "As":
        return AtomType(32)

    elif atom == "Se":
        return AtomType(33)

    elif atom == "Br":
        return AtomType(34)

    elif atom == "Kr":
        return AtomType(35)

    elif atom == "Rb":
        return AtomType(36)

    elif atom == "Sr":
        return AtomType(37)

    elif atom == "Zr":
        return AtomType(39)

    elif atom == "Nb":
        return AtomType(40)

    elif atom == "Mo":
        return AtomType(41)

    elif atom == "Tc":
        return AtomType(42)

    elif atom == "Ru":
        return AtomType(43)

    elif atom == "Rh":
        return AtomType(44)

    elif atom == "Pd":
        return AtomType(45)

    elif atom == "Ag":
        return AtomType(46)

    elif atom == "Cd":
        return AtomType(47)

    elif atom == "In":
        return AtomType(48)

    elif atom == "Sn":
        return AtomType(49)

    elif atom == "Sb":
        return AtomType(50)

    elif atom == "Te":
        return AtomType(51)

    elif atom == "I":
        return AtomType(52)

    elif atom == "Xe":
        return AtomType(53)

    elif atom == "Cs":
        return AtomType(54)

    elif atom == "Ba":
        return AtomType(55)

    elif atom == "La":
        return AtomType(56)

    elif atom == "Ce":
        return AtomType(57)

    elif atom == "Pr":
        return AtomType(58)

    elif atom == "Nd":
        return AtomType(59)

    elif atom == "Pm":
        return AtomType(60)

    elif atom == "Sm":
        return AtomType(61)

    elif atom == "Eu":
        return AtomType(62)

    elif atom == "Gd":
        return AtomType(63)

    elif atom == "Tb":
        return AtomType(64)

    elif atom == "Dy":
        return AtomType(65)

    elif atom == "Ho":
        return AtomType(66)

    elif atom == "Er":
        return AtomType(67)

    elif atom == "Tm":
        return AtomType(68)

    elif atom == "Yb":
        return AtomType(69)

    elif atom == "Lu":
        return AtomType(70)

    elif atom == "Hf":
        return AtomType(71)

    elif atom == "Ta":
        return AtomType(72)

    elif atom == "W":
        return AtomType(73)

    elif atom == "Re":
        return AtomType(74)

    elif atom == "Os":
        return AtomType(75)

    elif atom == "Ir":
        return AtomType(76)

    elif atom == "Pt":
        return AtomType(77)

    elif atom == "Au":
        return AtomType(78)

    elif atom == "Hg":
        return AtomType(79)

    elif atom == "Tl":
        return AtomType(80)

    elif atom == "Pb":
        return AtomType(81)

    elif atom == "Bi":
        return AtomType(82)

    elif atom == "Po":
        return AtomType(83)

    elif atom == "At":
        return AtomType(84)

    elif atom == "Rn":
        return AtomType(85)

    elif atom == "Fr":
        return AtomType(86)

    elif atom == "Ra":
        return AtomType(87)

    elif atom == "Ac":
        return AtomType(88)

    elif atom == "Th":
        return AtomType(89)

    elif atom == "Pa":
        return AtomType(90)

    elif atom == "U":
        return AtomType(91)

    elif atom == "Np":
        return AtomType(92)

    elif atom == "Pu":
        return AtomType(93)

    elif atom == "Am":
        return AtomType(94)

    elif atom == "Cm":
        return AtomType(95)

    elif atom == "Bk":
        return AtomType(96)

    elif atom == "Cf":
        return AtomType(97)

    elif atom == "Es":
        return AtomType(98)

    elif atom == "Fm":
        return AtomType(99)

    elif atom == "Md":
        return AtomType(100)

    elif atom == "No":
        return AtomType(101)

    elif atom == "Lr":
        return AtomType(102)

    elif atom == "Rf":
        return AtomType(103)

    elif atom == "Db":
        return AtomType(104)

    elif atom == "Sg":
        return AtomType(105)

    elif atom == "Bh":
        return AtomType(106)

    elif atom == "Hs":
        return AtomType(107)

    elif atom == "Mt":
        return AtomType(108)

    elif atom == "Ds":
        return AtomType(109)

    elif atom == "Rg":
        return AtomType(110)

    elif atom == "Cn":
        return AtomType(111)

    elif atom == "Nh":
        return AtomType(112)

    elif atom == "Fl":
        return AtomType(113)

    elif atom == "Mc":
        return AtomType(114)

    elif atom == "Lv":
        return AtomType(115)

    elif atom == "Ts":
        return AtomType(116)

    elif atom == "Og":
        return AtomType(117)

    elif atom == "Y":
        return AtomType(38)

    else: 
        return AtomType(118)



def try_cast_to_bond(bond: str) -> BondType:
    if bond == "1":
        return BondType(0)

    elif bond == "2":
        return BondType(1)

    elif bond == "3":
        return BondType(2)

    elif bond == "4":
        return BondType(3)

    else: 
        return BondType(4)



def try_cast_to_float(s_1: str) -> float:
    try: 
        return parse(s_1)

    except Exception as match_value:
        raise match_value



def try_cast_to_int(s_1: str) -> int:
    try: 
        return parse_1(s_1, 511, False, 32)

    except Exception as match_value:
        raise match_value



def parse_sdf(sdf: str) -> Array[Molecule]:
    atoms: FSharpList[AtomInfo] = empty()
    atom_count: int = 0
    bonds: FSharpList[BondInfo] = empty()
    bond_count: int = 0
    def _arrow55(__unit: None=None, sdf: str=sdf) -> IEnumerable_1[Molecule]:
        def _arrow54(line: str) -> IEnumerable_1[Molecule]:
            nonlocal atom_count, atoms, bond_count, bonds
            match_value: str = line
            (pattern_matching_result, symbol, x, y, z) = (None, None, None, None, None)
            if (match_value.find("$$$$") >= 0) == True:
                pattern_matching_result = 0

            else: 
                active_pattern_result: Optional[FSharpList[str]] = _007CAtomLine_007C__007C(match_value)
                if active_pattern_result is not None:
                    if not is_empty(active_pattern_result):
                        if not is_empty(tail(active_pattern_result)):
                            if not is_empty(tail(tail(active_pattern_result))):
                                if not is_empty(tail(tail(tail(active_pattern_result)))):
                                    if is_empty(tail(tail(tail(tail(active_pattern_result))))):
                                        pattern_matching_result = 1
                                        symbol = head(tail(tail(tail(active_pattern_result))))
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
                def _arrow53(__unit: None=None) -> IEnumerable_1[Molecule]:
                    nonlocal atom_count, atoms
                    atom_count = 0
                    atoms = empty()
                    return empty_1()

                return append_1(singleton_1(Molecule(atoms, bonds)), delay(_arrow53))

            elif pattern_matching_result == 1:
                atom_count = (atom_count + 1) or 0
                atom_type: AtomType = try_cast_to_atom(symbol)
                center: Point3D = Point3D(try_cast_to_float(x), try_cast_to_float(y), try_cast_to_float(z))
                radius_1: float = normalize_radius(AtomGeomStyle(0), get_atom_radius(AtomGeomStyle(0), atom_type))
                atom: AtomInfo = create_atom(atom_count, atom_type, center, radius_1)
                atoms = append(atoms, singleton(atom))
                return empty_1()

            elif pattern_matching_result == 2:
                (pattern_matching_result_1, bond_type, e, s_1) = (None, None, None, None)
                active_pattern_result_1: Optional[FSharpList[str]] = _007CBondLine_007C__007C(match_value)
                if active_pattern_result_1 is not None:
                    if not is_empty(active_pattern_result_1):
                        if not is_empty(tail(active_pattern_result_1)):
                            if not is_empty(tail(tail(active_pattern_result_1))):
                                if is_empty(tail(tail(tail(active_pattern_result_1)))):
                                    pattern_matching_result_1 = 0
                                    bond_type = head(tail(tail(active_pattern_result_1)))
                                    e = head(tail(active_pattern_result_1))
                                    s_1 = head(active_pattern_result_1)

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
                    bond_type_1: BondType = try_cast_to_bond(bond_type)
                    s_2: int = try_cast_to_int(s_1) or 0
                    e_1: int = try_cast_to_int(e) or 0
                    bond_count = (bond_count + 1) or 0
                    bond: BondInfo = BondInfo(bond_count, s_2, e_1, bond_type_1, 1.0)
                    bond_count = (bond_count + 1) or 0
                    rev_bond: BondInfo = BondInfo(bond_count, e_1, s_2, bond_type_1, 1.0)
                    bonds = append(bonds, of_array([bond, rev_bond]))
                    return empty_1()

                elif pattern_matching_result_1 == 1:
                    return empty_1()



        return collect(_arrow54, sdf.split("\n"))

    return to_array(delay(_arrow55))


__all__ = ["ParserError_reflection", "ParserError__Equals_229D3F39", "s", "d", "d_cap", "w_cap", "atom_line", "_007CAtomLine_007C__007C", "bond_line", "_007CBondLine_007C__007C", "try_cast_to_atom", "try_cast_to_bond", "try_cast_to_float", "try_cast_to_int", "parse_sdf"]

