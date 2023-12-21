from math import pow
from .fable_modules.fable_library.double import sqrt
from .fable_modules.fable_library.util import round as round_1

def round(digits: int, v: float) -> float:
    return round_1(v, digits)


def abs_1(v: float) -> float:
    return sqrt(pow(v, 2.0))


def clamp(lower_bound: float, upper_bound: float, v: float) -> float:
    if v < lower_bound:
        return lower_bound

    elif v > upper_bound:
        return upper_bound

    else: 
        return v



__all__ = ["round", "abs_1", "clamp"]

