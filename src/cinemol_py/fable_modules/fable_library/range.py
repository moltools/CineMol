from __future__ import annotations

from collections.abc import Callable
from decimal import Decimal
from typing import Any, TypeVar

from .big_int import from_zero, op_addition
from .char import char_code_at
from .decimal_ import from_parts
from .decimal_ import op_addition as op_addition_1
from .long import op_addition as op_addition_2
from .seq import delay, unfold
from .types import int64, uint8, uint64
from .util import IEnumerable_1, compare


_T = TypeVar("_T")


def make_range_step_function(
    step: _T, stop: _T, zero: _T, add: Callable[[_T, _T], _T]
) -> Callable[[_T], tuple[_T, _T] | None]:
    step_compared_with_zero: int = compare(step, zero) or 0
    if step_compared_with_zero == 0:
        raise Exception("The step of a range cannot be zero")

    step_greater_than_zero: bool = step_compared_with_zero > 0

    def _arrow128(
        x: _T | None = None, step: Any = step, stop: Any = stop, zero: Any = zero, add: Any = add
    ) -> tuple[_T, _T] | None:
        compared_with_last: int = compare(x, stop) or 0
        return (
            ((x, add(x, step)))
            if (
                True
                if ((compared_with_last <= 0) if step_greater_than_zero else False)
                else ((compared_with_last >= 0) if (not step_greater_than_zero) else False)
            )
            else None
        )

    return _arrow128


def integral_range_step(start: _T, step: _T, stop: _T, zero: _T, add: Callable[[_T, _T], _T]) -> IEnumerable_1[_T]:
    step_fn: Callable[[_T], tuple[_T, _T] | None] = make_range_step_function(step, stop, zero, add)

    def _arrow129(
        __unit: None = None, start: Any = start, step: Any = step, stop: Any = stop, zero: Any = zero, add: Any = add
    ) -> IEnumerable_1[_T]:
        return unfold(step_fn, start)

    return delay(_arrow129)


def range_big_int(start: int, step: int, stop: int) -> IEnumerable_1[int]:
    def _arrow130(x: int, y: int, start: Any = start, step: Any = step, stop: Any = stop) -> int:
        return op_addition(x, y)

    return integral_range_step(start, step, stop, from_zero(), _arrow130)


def range_decimal(start: Decimal, step: Decimal, stop: Decimal) -> IEnumerable_1[Decimal]:
    def _arrow131(x: Decimal, y: Decimal, start: Any = start, step: Any = step, stop: Any = stop) -> Decimal:
        return op_addition_1(x, y)

    return integral_range_step(start, step, stop, from_parts(0, 0, 0, False, uint8(0)), _arrow131)


def range_double(start: float, step: float, stop: float) -> IEnumerable_1[float]:
    def _arrow132(x: float, y: float, start: Any = start, step: Any = step, stop: Any = stop) -> float:
        return x + y

    return integral_range_step(start, step, stop, 0.0, _arrow132)


def range_int64(start: int64, step: int64, stop: int64) -> IEnumerable_1[int64]:
    def _arrow133(x: int64, y: int64, start: Any = start, step: Any = step, stop: Any = stop) -> int64:
        return op_addition_2(x, y)

    return integral_range_step(start, step, stop, int64(0), _arrow133)


def range_uint64(start: uint64, step: uint64, stop: uint64) -> IEnumerable_1[uint64]:
    def _arrow134(x: uint64, y: uint64, start: Any = start, step: Any = step, stop: Any = stop) -> uint64:
        return op_addition_2(x, y)

    return integral_range_step(start, step, stop, uint64(0), _arrow134)


def range_char(start: str, stop: str) -> IEnumerable_1[str]:
    int_stop: int = char_code_at(stop, 0) or 0

    def _arrow135(__unit: None = None, start: Any = start, stop: Any = stop) -> IEnumerable_1[str]:
        def step_fn(c: int) -> tuple[str, int] | None:
            if c <= int_stop:
                return (chr(c), c + 1)

            else:
                return None

        return unfold(step_fn, char_code_at(start, 0))

    return delay(_arrow135)


__all__ = [
    "make_range_step_function",
    "integral_range_step",
    "range_big_int",
    "range_decimal",
    "range_double",
    "range_int64",
    "range_uint64",
    "range_char",
]
