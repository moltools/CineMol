from decimal import Decimal
from typing import (TypeVar, Callable, Any, Optional, Tuple)
from .big_int import (from_zero, op_addition)
from .char import char_code_at
from .decimal import (from_parts, op_addition as op_addition_1)
from .long import op_addition as op_addition_2
from .seq import (delay, unfold)
from .types import (uint8, int64, uint64)
from .util import (compare, IEnumerable_1)

_T = TypeVar("_T")

def make_range_step_function(step: _T, stop: _T, zero: _T, add: Callable[[_T, _T], _T]) -> Callable[[_T], Optional[Tuple[_T, _T]]]:
    step_compared_with_zero: int = compare(step, zero) or 0
    if step_compared_with_zero == 0:
        raise Exception("The step of a range cannot be zero")

    step_greater_than_zero: bool = step_compared_with_zero > 0
    def _arrow125(x: Optional[_T]=None, step: _T=step, stop: _T=stop, zero: _T=zero, add: Any=add) -> Optional[Tuple[_T, _T]]:
        compared_with_last: int = compare(x, stop) or 0
        return ((x, add(x, step))) if (True if ((compared_with_last <= 0) if step_greater_than_zero else False) else ((compared_with_last >= 0) if (not step_greater_than_zero) else False)) else None

    return _arrow125


def integral_range_step(start: _T, step: _T, stop: _T, zero: _T, add: Callable[[_T, _T], _T]) -> IEnumerable_1[_T]:
    step_fn: Callable[[_T], Optional[Tuple[_T, _T]]] = make_range_step_function(step, stop, zero, add)
    def _arrow126(__unit: None=None, start: _T=start, step: _T=step, stop: _T=stop, zero: _T=zero, add: Any=add) -> IEnumerable_1[_T]:
        return unfold(step_fn, start)

    return delay(_arrow126)


def range_big_int(start: int, step: int, stop: int) -> IEnumerable_1[int]:
    def _arrow127(x: int, y: int, start: int=start, step: int=step, stop: int=stop) -> int:
        return op_addition(x, y)

    return integral_range_step(start, step, stop, from_zero(), _arrow127)


def range_decimal(start: Decimal, step: Decimal, stop: Decimal) -> IEnumerable_1[Decimal]:
    def _arrow128(x: Decimal, y: Decimal, start: Decimal=start, step: Decimal=step, stop: Decimal=stop) -> Decimal:
        return op_addition_1(x, y)

    return integral_range_step(start, step, stop, from_parts(0, 0, 0, False, uint8(0)), _arrow128)


def range_double(start: float, step: float, stop: float) -> IEnumerable_1[float]:
    def _arrow129(x: float, y: float, start: float=start, step: float=step, stop: float=stop) -> float:
        return x + y

    return integral_range_step(start, step, stop, 0.0, _arrow129)


def range_int64(start: int64, step: int64, stop: int64) -> IEnumerable_1[int64]:
    def _arrow130(x: int64, y: int64, start: int64=start, step: int64=step, stop: int64=stop) -> int64:
        return op_addition_2(x, y)

    return integral_range_step(start, step, stop, int64(0), _arrow130)


def range_uint64(start: uint64, step: uint64, stop: uint64) -> IEnumerable_1[uint64]:
    def _arrow131(x: uint64, y: uint64, start: uint64=start, step: uint64=step, stop: uint64=stop) -> uint64:
        return op_addition_2(x, y)

    return integral_range_step(start, step, stop, uint64(0), _arrow131)


def range_char(start: str, stop: str) -> IEnumerable_1[str]:
    int_stop: int = char_code_at(stop, 0) or 0
    def step_fn(c: int, start: str=start, stop: str=stop) -> Optional[Tuple[str, int]]:
        if c <= int_stop:
            return (chr(c), c + 1)

        else: 
            return None


    def _arrow132(__unit: None=None, start: str=start, stop: str=stop) -> IEnumerable_1[str]:
        return unfold(step_fn, char_code_at(start, 0))

    return delay(_arrow132)


__all__ = ["make_range_step_function", "integral_range_step", "range_big_int", "range_decimal", "range_double", "range_int64", "range_uint64", "range_char"]

