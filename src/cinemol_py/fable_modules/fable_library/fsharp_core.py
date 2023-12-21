from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

from .fsharp_collections import ComparisonIdentity_Structural, HashIdentity_Structural
from .system_text import StringBuilder, StringBuilder__Append_Z721C83C5
from .util import IComparer_1, IEqualityComparer, IEqualityComparer_1, dispose, equals, ignore, structural_hash


_T = TypeVar("_T")

_R = TypeVar("_R")

__B = TypeVar("__B")

__A = TypeVar("__A")


class ObjectExpr1(IEqualityComparer):
    def System_Collections_IEqualityComparer_Equals541DA560(self, x: Any = None, y: Any = None) -> bool:
        return equals(x, y)

    def System_Collections_IEqualityComparer_GetHashCode4E60E31B(self, x_1: Any = None) -> int:
        return structural_hash(x_1)


LanguagePrimitives_GenericEqualityComparer: IEqualityComparer = ObjectExpr1()


class ObjectExpr2(IEqualityComparer):
    def System_Collections_IEqualityComparer_Equals541DA560(self, x: Any = None, y: Any = None) -> bool:
        return equals(x, y)

    def System_Collections_IEqualityComparer_GetHashCode4E60E31B(self, x_1: Any = None) -> int:
        return structural_hash(x_1)


LanguagePrimitives_GenericEqualityERComparer: IEqualityComparer = ObjectExpr2()


def LanguagePrimitives_FastGenericComparer(__unit: None = None) -> IComparer_1[Any]:
    return ComparisonIdentity_Structural()


def LanguagePrimitives_FastGenericComparerFromTable(__unit: None = None) -> IComparer_1[Any]:
    return ComparisonIdentity_Structural()


def LanguagePrimitives_FastGenericEqualityComparer(__unit: None = None) -> IEqualityComparer_1[Any]:
    return HashIdentity_Structural()


def LanguagePrimitives_FastGenericEqualityComparerFromTable(__unit: None = None) -> IEqualityComparer_1[Any]:
    return HashIdentity_Structural()


def Operators_Failure(message: str) -> Exception:
    return Exception(message)


def Operators_FailurePattern(exn: Exception) -> str | None:
    return str(exn)


def Operators_NullArg(x: str) -> Any:
    raise Exception(x)


def Operators_Using(resource: _T, action: Callable[[_T], _R]) -> _R:
    try:
        return action(resource)

    finally:
        if equals(resource, None):
            pass

        else:
            copy_of_struct: _T = resource
            dispose(copy_of_struct)


def Operators_Lock(_lockObj: Any, action: Callable[[], __B]) -> __B:
    return action(None)


def ExtraTopLevelOperators_LazyPattern(input: Any) -> __A:
    return input.Value


def PrintfModule_PrintFormatToStringBuilderThen(
    continuation: Callable[[], __A], builder: StringBuilder, format: Any
) -> __B:
    def append(s: str, continuation: Any = continuation, builder: Any = builder, format: Any = format) -> __A:
        ignore(StringBuilder__Append_Z721C83C5(builder, s))
        return continuation(None)

    return format.cont(append)


def PrintfModule_PrintFormatToStringBuilder(builder: StringBuilder, format: Any) -> __A:
    def _arrow5(__unit: None = None, builder: Any = builder, format: Any = format) -> None:
        ignore(None)

    return PrintfModule_PrintFormatToStringBuilderThen(_arrow5, builder, format)


__all__ = [
    "LanguagePrimitives_GenericEqualityComparer",
    "LanguagePrimitives_GenericEqualityERComparer",
    "LanguagePrimitives_FastGenericComparer",
    "LanguagePrimitives_FastGenericComparerFromTable",
    "LanguagePrimitives_FastGenericEqualityComparer",
    "LanguagePrimitives_FastGenericEqualityComparerFromTable",
    "Operators_Failure",
    "Operators_FailurePattern",
    "Operators_NullArg",
    "Operators_Using",
    "Operators_Lock",
    "ExtraTopLevelOperators_LazyPattern",
    "PrintfModule_PrintFormatToStringBuilderThen",
    "PrintfModule_PrintFormatToStringBuilder",
]
