from __future__ import annotations
from typing import Any
from .reflection import (TypeInfo, class_type)
from .string import (join, is_null_or_empty, substring, format)
from .types import (Array, to_string)
from .util import (int32_to_string, clear)

def _expr0() -> TypeInfo:
    return class_type("System.Text.StringBuilder", None, StringBuilder)


class StringBuilder:
    def __init__(self, value: str, capacity: int) -> None:
        self.buf: Array[str] = []
        if not is_null_or_empty(value):
            (self.buf.append(value))


    def __str__(self, __unit: None=None) -> str:
        __: StringBuilder = self
        return join("", __.buf)


StringBuilder_reflection = _expr0

def StringBuilder__ctor_Z18115A39(value: str, capacity: int) -> StringBuilder:
    return StringBuilder(value, capacity)


def StringBuilder__ctor_Z524259A4(capacity: int) -> StringBuilder:
    return StringBuilder__ctor_Z18115A39("", capacity)


def StringBuilder__ctor_Z721C83C5(value: str) -> StringBuilder:
    return StringBuilder__ctor_Z18115A39(value, 16)


def StringBuilder__ctor(__unit: None=None) -> StringBuilder:
    return StringBuilder__ctor_Z18115A39("", 16)


def StringBuilder__Append_Z721C83C5(x: StringBuilder, s: str) -> StringBuilder:
    (x.buf.append(s))
    return x


def StringBuilder__Append_487EF8FB(x: StringBuilder, s: str, start_index: int, count: int) -> StringBuilder:
    (x.buf.append(substring(s, start_index, count)))
    return x


def StringBuilder__Append_244C7CD6(x: StringBuilder, c: str) -> StringBuilder:
    (x.buf.append(c))
    return x


def StringBuilder__Append_Z524259A4(x: StringBuilder, o: int) -> StringBuilder:
    (x.buf.append(int32_to_string(o)))
    return x


def StringBuilder__Append_5E38073B(x: StringBuilder, o: float) -> StringBuilder:
    (x.buf.append(to_string(o)))
    return x


def StringBuilder__Append_Z1FBCCD16(x: StringBuilder, o: bool) -> StringBuilder:
    (x.buf.append(to_string(o)))
    return x


def StringBuilder__Append_4E60E31B(x: StringBuilder, o: Any=None) -> StringBuilder:
    (x.buf.append(to_string(o)))
    return x


def StringBuilder__Append_Z372E4D23(x: StringBuilder, cs: Array[str]) -> StringBuilder:
    (x.buf.append(''.join(cs)))
    return x


def StringBuilder__Append_43A65C09(x: StringBuilder, s: StringBuilder) -> StringBuilder:
    (x.buf.append(to_string(s)))
    return x


def StringBuilder__AppendFormat_433E080(x: StringBuilder, fmt: str, o: Any=None) -> StringBuilder:
    (x.buf.append(format(fmt, o)))
    return x


def StringBuilder__AppendFormat_Z696D8D1B(x: StringBuilder, provider: IFormatProvider, fmt: str, o: Any=None) -> StringBuilder:
    (x.buf.append(format(provider, fmt, o)))
    return x


def StringBuilder__AppendLine(x: StringBuilder) -> StringBuilder:
    (x.buf.append("\n"))
    return x


def StringBuilder__AppendLine_Z721C83C5(x: StringBuilder, s: str) -> StringBuilder:
    (x.buf.append(s))
    (x.buf.append("\n"))
    return x


def StringBuilder__get_Length(x: StringBuilder) -> int:
    len_1: int = 0
    for i in range(len(x.buf) - 1, 0 - 1, -1):
        len_1 = (len_1 + len(x.buf[i])) or 0
    return len_1


def StringBuilder__ToString_Z37302880(x: StringBuilder, first_index: int, length: int) -> str:
    return substring(to_string(x), first_index, length)


def StringBuilder__Clear(x: StringBuilder) -> StringBuilder:
    clear(x.buf)
    return x


__all__ = ["StringBuilder_reflection", "StringBuilder__ctor_Z524259A4", "StringBuilder__ctor_Z721C83C5", "StringBuilder__ctor", "StringBuilder__Append_Z721C83C5", "StringBuilder__Append_487EF8FB", "StringBuilder__Append_244C7CD6", "StringBuilder__Append_Z524259A4", "StringBuilder__Append_5E38073B", "StringBuilder__Append_Z1FBCCD16", "StringBuilder__Append_4E60E31B", "StringBuilder__Append_Z372E4D23", "StringBuilder__Append_43A65C09", "StringBuilder__AppendFormat_433E080", "StringBuilder__AppendFormat_Z696D8D1B", "StringBuilder__AppendLine", "StringBuilder__AppendLine_Z721C83C5", "StringBuilder__get_Length", "StringBuilder__ToString_Z37302880", "StringBuilder__Clear"]

