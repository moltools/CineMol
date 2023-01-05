from math import pow
from .fable_modules.fable_library.types import to_string
from .fable_modules.fable_library.util import (round as round_1, int32_to_string)

def round(digits: int, f: float) -> float:
    return round_1(f, digits)


def float_to_str(f: float) -> str:
    return to_string(f)


def int_to_str(d: int) -> str:
    return str(d)


def abs_1(f: float) -> float:
    return pow(pow(f, 2.0), 0.5)


__all__ = ["round", "float_to_str", "int_to_str", "abs_1"]

