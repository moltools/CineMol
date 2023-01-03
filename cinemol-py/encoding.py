from __future__ import annotations
from typing import (Any, List)
from .fable_modules.fable_library.array import map
from .fable_modules.fable_library.reflection import (TypeInfo, class_type, union_type)
from .fable_modules.fable_library.seq import to_array
from .fable_modules.fable_library.string import to_base64string as to_base64string_1
from .fable_modules.fable_library.types import (FSharpException, Array, Union, uint8, Uint8Array)
from .fable_modules.fable_library.util import equals

def _expr1() -> TypeInfo:
    return class_type("Client.CineMol.Encoding.UnexpectedCharToEncode", None, UnexpectedCharToEncode, class_type("System.Exception"))


class UnexpectedCharToEncode(FSharpException):
    def __init__(self, Data0: str) -> None:
        super().__init__()
        self.Data0 = Data0


UnexpectedCharToEncode_reflection = _expr1

def UnexpectedCharToEncode__Equals_229D3F39(this: Exception, obj: Exception) -> bool:
    if not equals(this, None):
        if not equals(obj, None):
            if isinstance(obj, UnexpectedCharToEncode):
                return this.Data0 == obj.Data0

            else: 
                return False


        else: 
            return False


    elif not equals(obj, None):
        return False

    else: 
        return True



def _expr2() -> TypeInfo:
    return union_type("Client.CineMol.Encoding.Encoding", [], Encoding, lambda: [[]])


class Encoding(Union):
    def __init__(self, tag: int, *fields: Any) -> None:
        super().__init__()
        self.tag: int = tag or 0
        self.fields: Array[Any] = list(fields)

    @staticmethod
    def cases() -> List[str]:
        return ["ISO_8859_1"]


Encoding_reflection = _expr2

def Encoding__GetEncoding_244C7CD6(x: Encoding, c: str) -> uint8:
    (pattern_matching_result,) = (None,)
    if c == "\t":
        pattern_matching_result = 0

    elif c == "\n":
        pattern_matching_result = 0

    elif c == " ":
        pattern_matching_result = 0

    elif c == "!":
        pattern_matching_result = 1

    elif c == "\"":
        pattern_matching_result = 2

    elif c == "#":
        pattern_matching_result = 3

    elif c == "$":
        pattern_matching_result = 4

    elif c == "%":
        pattern_matching_result = 5

    elif c == "&":
        pattern_matching_result = 6

    elif c == "\'":
        pattern_matching_result = 7

    elif c == "(":
        pattern_matching_result = 8

    elif c == ")":
        pattern_matching_result = 9

    elif c == "*":
        pattern_matching_result = 10

    elif c == "+":
        pattern_matching_result = 11

    elif c == ",":
        pattern_matching_result = 12

    elif c == "-":
        pattern_matching_result = 13

    elif c == ".":
        pattern_matching_result = 14

    elif c == "/":
        pattern_matching_result = 15

    elif c == "0":
        pattern_matching_result = 16

    elif c == "1":
        pattern_matching_result = 17

    elif c == "2":
        pattern_matching_result = 18

    elif c == "3":
        pattern_matching_result = 19

    elif c == "4":
        pattern_matching_result = 20

    elif c == "5":
        pattern_matching_result = 21

    elif c == "6":
        pattern_matching_result = 22

    elif c == "7":
        pattern_matching_result = 23

    elif c == "8":
        pattern_matching_result = 24

    elif c == "9":
        pattern_matching_result = 25

    elif c == ":":
        pattern_matching_result = 26

    elif c == ";":
        pattern_matching_result = 27

    elif c == "<":
        pattern_matching_result = 28

    elif c == "=":
        pattern_matching_result = 29

    elif c == ">":
        pattern_matching_result = 30

    elif c == "?":
        pattern_matching_result = 31

    elif c == "@":
        pattern_matching_result = 32

    elif c == "A":
        pattern_matching_result = 33

    elif c == "B":
        pattern_matching_result = 34

    elif c == "C":
        pattern_matching_result = 35

    elif c == "D":
        pattern_matching_result = 36

    elif c == "E":
        pattern_matching_result = 37

    elif c == "F":
        pattern_matching_result = 38

    elif c == "G":
        pattern_matching_result = 39

    elif c == "H":
        pattern_matching_result = 40

    elif c == "I":
        pattern_matching_result = 41

    elif c == "J":
        pattern_matching_result = 42

    elif c == "K":
        pattern_matching_result = 43

    elif c == "L":
        pattern_matching_result = 44

    elif c == "M":
        pattern_matching_result = 45

    elif c == "N":
        pattern_matching_result = 46

    elif c == "O":
        pattern_matching_result = 47

    elif c == "P":
        pattern_matching_result = 48

    elif c == "Q":
        pattern_matching_result = 49

    elif c == "R":
        pattern_matching_result = 50

    elif c == "S":
        pattern_matching_result = 51

    elif c == "T":
        pattern_matching_result = 52

    elif c == "U":
        pattern_matching_result = 53

    elif c == "V":
        pattern_matching_result = 54

    elif c == "W":
        pattern_matching_result = 55

    elif c == "X":
        pattern_matching_result = 56

    elif c == "Y":
        pattern_matching_result = 57

    elif c == "Z":
        pattern_matching_result = 58

    elif c == "[":
        pattern_matching_result = 59

    elif c == "\\":
        pattern_matching_result = 60

    elif c == "]":
        pattern_matching_result = 61

    elif c == "^":
        pattern_matching_result = 62

    elif c == "_":
        pattern_matching_result = 63

    elif c == "`":
        pattern_matching_result = 64

    elif c == "a":
        pattern_matching_result = 65

    elif c == "b":
        pattern_matching_result = 66

    elif c == "c":
        pattern_matching_result = 67

    elif c == "d":
        pattern_matching_result = 68

    elif c == "e":
        pattern_matching_result = 69

    elif c == "f":
        pattern_matching_result = 70

    elif c == "g":
        pattern_matching_result = 71

    elif c == "h":
        pattern_matching_result = 72

    elif c == "i":
        pattern_matching_result = 73

    elif c == "j":
        pattern_matching_result = 74

    elif c == "k":
        pattern_matching_result = 75

    elif c == "l":
        pattern_matching_result = 76

    elif c == "m":
        pattern_matching_result = 77

    elif c == "n":
        pattern_matching_result = 78

    elif c == "o":
        pattern_matching_result = 79

    elif c == "p":
        pattern_matching_result = 80

    elif c == "q":
        pattern_matching_result = 81

    elif c == "r":
        pattern_matching_result = 82

    elif c == "s":
        pattern_matching_result = 83

    elif c == "t":
        pattern_matching_result = 84

    elif c == "u":
        pattern_matching_result = 85

    elif c == "v":
        pattern_matching_result = 86

    elif c == "w":
        pattern_matching_result = 87

    elif c == "x":
        pattern_matching_result = 88

    elif c == "y":
        pattern_matching_result = 89

    elif c == "z":
        pattern_matching_result = 90

    elif c == "{":
        pattern_matching_result = 91

    elif c == "|":
        pattern_matching_result = 92

    elif c == "}":
        pattern_matching_result = 93

    elif c == "~":
        pattern_matching_result = 94

    else: 
        pattern_matching_result = 95

    if pattern_matching_result == 0:
        return uint8(32)

    elif pattern_matching_result == 1:
        return uint8(33)

    elif pattern_matching_result == 2:
        return uint8(34)

    elif pattern_matching_result == 3:
        return uint8(35)

    elif pattern_matching_result == 4:
        return uint8(36)

    elif pattern_matching_result == 5:
        return uint8(37)

    elif pattern_matching_result == 6:
        return uint8(38)

    elif pattern_matching_result == 7:
        return uint8(39)

    elif pattern_matching_result == 8:
        return uint8(40)

    elif pattern_matching_result == 9:
        return uint8(41)

    elif pattern_matching_result == 10:
        return uint8(42)

    elif pattern_matching_result == 11:
        return uint8(43)

    elif pattern_matching_result == 12:
        return uint8(44)

    elif pattern_matching_result == 13:
        return uint8(45)

    elif pattern_matching_result == 14:
        return uint8(46)

    elif pattern_matching_result == 15:
        return uint8(47)

    elif pattern_matching_result == 16:
        return uint8(48)

    elif pattern_matching_result == 17:
        return uint8(49)

    elif pattern_matching_result == 18:
        return uint8(50)

    elif pattern_matching_result == 19:
        return uint8(51)

    elif pattern_matching_result == 20:
        return uint8(52)

    elif pattern_matching_result == 21:
        return uint8(53)

    elif pattern_matching_result == 22:
        return uint8(54)

    elif pattern_matching_result == 23:
        return uint8(55)

    elif pattern_matching_result == 24:
        return uint8(56)

    elif pattern_matching_result == 25:
        return uint8(57)

    elif pattern_matching_result == 26:
        return uint8(58)

    elif pattern_matching_result == 27:
        return uint8(59)

    elif pattern_matching_result == 28:
        return uint8(60)

    elif pattern_matching_result == 29:
        return uint8(61)

    elif pattern_matching_result == 30:
        return uint8(62)

    elif pattern_matching_result == 31:
        return uint8(63)

    elif pattern_matching_result == 32:
        return uint8(64)

    elif pattern_matching_result == 33:
        return uint8(65)

    elif pattern_matching_result == 34:
        return uint8(66)

    elif pattern_matching_result == 35:
        return uint8(67)

    elif pattern_matching_result == 36:
        return uint8(68)

    elif pattern_matching_result == 37:
        return uint8(69)

    elif pattern_matching_result == 38:
        return uint8(70)

    elif pattern_matching_result == 39:
        return uint8(71)

    elif pattern_matching_result == 40:
        return uint8(72)

    elif pattern_matching_result == 41:
        return uint8(73)

    elif pattern_matching_result == 42:
        return uint8(74)

    elif pattern_matching_result == 43:
        return uint8(75)

    elif pattern_matching_result == 44:
        return uint8(76)

    elif pattern_matching_result == 45:
        return uint8(77)

    elif pattern_matching_result == 46:
        return uint8(78)

    elif pattern_matching_result == 47:
        return uint8(79)

    elif pattern_matching_result == 48:
        return uint8(80)

    elif pattern_matching_result == 49:
        return uint8(81)

    elif pattern_matching_result == 50:
        return uint8(82)

    elif pattern_matching_result == 51:
        return uint8(83)

    elif pattern_matching_result == 52:
        return uint8(84)

    elif pattern_matching_result == 53:
        return uint8(85)

    elif pattern_matching_result == 54:
        return uint8(86)

    elif pattern_matching_result == 55:
        return uint8(87)

    elif pattern_matching_result == 56:
        return uint8(88)

    elif pattern_matching_result == 57:
        return uint8(89)

    elif pattern_matching_result == 58:
        return uint8(90)

    elif pattern_matching_result == 59:
        return uint8(91)

    elif pattern_matching_result == 60:
        return uint8(92)

    elif pattern_matching_result == 61:
        return uint8(93)

    elif pattern_matching_result == 62:
        return uint8(94)

    elif pattern_matching_result == 63:
        return uint8(95)

    elif pattern_matching_result == 64:
        return uint8(96)

    elif pattern_matching_result == 65:
        return uint8(97)

    elif pattern_matching_result == 66:
        return uint8(98)

    elif pattern_matching_result == 67:
        return uint8(99)

    elif pattern_matching_result == 68:
        return uint8(100)

    elif pattern_matching_result == 69:
        return uint8(101)

    elif pattern_matching_result == 70:
        return uint8(102)

    elif pattern_matching_result == 71:
        return uint8(103)

    elif pattern_matching_result == 72:
        return uint8(104)

    elif pattern_matching_result == 73:
        return uint8(105)

    elif pattern_matching_result == 74:
        return uint8(106)

    elif pattern_matching_result == 75:
        return uint8(107)

    elif pattern_matching_result == 76:
        return uint8(108)

    elif pattern_matching_result == 77:
        return uint8(109)

    elif pattern_matching_result == 78:
        return uint8(110)

    elif pattern_matching_result == 79:
        return uint8(111)

    elif pattern_matching_result == 80:
        return uint8(112)

    elif pattern_matching_result == 81:
        return uint8(113)

    elif pattern_matching_result == 82:
        return uint8(114)

    elif pattern_matching_result == 83:
        return uint8(115)

    elif pattern_matching_result == 84:
        return uint8(116)

    elif pattern_matching_result == 85:
        return uint8(117)

    elif pattern_matching_result == 86:
        return uint8(118)

    elif pattern_matching_result == 87:
        return uint8(119)

    elif pattern_matching_result == 88:
        return uint8(120)

    elif pattern_matching_result == 89:
        return uint8(121)

    elif pattern_matching_result == 90:
        return uint8(122)

    elif pattern_matching_result == 91:
        return uint8(123)

    elif pattern_matching_result == 92:
        return uint8(124)

    elif pattern_matching_result == 93:
        return uint8(125)

    elif pattern_matching_result == 94:
        return uint8(126)

    elif pattern_matching_result == 95:
        raise UnexpectedCharToEncode(((("Unknown character to encode " + str(x)) + ": ") + str(c)) + "")



def encode(encoding: Encoding, src: str) -> bytearray:
    def mapping(c: str, encoding: Encoding=encoding, src: str=src) -> uint8:
        return Encoding__GetEncoding_244C7CD6(encoding, c)

    return map(mapping, to_array(src), Uint8Array)


def to_base64string(to_encode: str) -> str:
    return to_base64string_1(encode(Encoding(0), to_encode))


__all__ = ["UnexpectedCharToEncode_reflection", "UnexpectedCharToEncode__Equals_229D3F39", "Encoding_reflection", "Encoding__GetEncoding_244C7CD6", "encode", "to_base64string"]

