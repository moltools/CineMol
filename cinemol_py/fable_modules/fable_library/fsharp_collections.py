from typing import (TypeVar, Callable, Any, Optional)
from .util import (IEqualityComparer_1, structural_hash, equals, physical_hash, IComparer_1, compare)

_T = TypeVar("_T")

_T_ = TypeVar("_T_")

def HashIdentity_FromFunctions(hash_1: Callable[[_T], int], eq: Callable[[_T, _T], bool]) -> IEqualityComparer_1[Any]:
    class ObjectExpr1(IEqualityComparer_1[Any]):
        def Equals(self, x: _T_, y: _T_, hash_1: Any=hash_1, eq: Any=eq) -> bool:
            return eq(x, y)

        def GetHashCode(self, x_1: Optional[_T_]=None, hash_1: Any=hash_1, eq: Any=eq) -> int:
            return hash_1(x_1)

    return ObjectExpr1()


def HashIdentity_Structural(__unit: None=None) -> IEqualityComparer_1[Any]:
    return HashIdentity_FromFunctions(structural_hash, equals)


def HashIdentity_Reference(__unit: None=None) -> IEqualityComparer_1[Any]:
    def _arrow2(e: _T, e_1: _T) -> bool:
        return e == e_1

    return HashIdentity_FromFunctions(physical_hash, _arrow2)


def ComparisonIdentity_FromFunction(comparer: Callable[[_T, _T], int]) -> IComparer_1[_T]:
    class ObjectExpr3(IComparer_1[_T_]):
        def Compare(self, x: _T_, y: _T_, comparer: Any=comparer) -> int:
            return comparer(x, y)

    return ObjectExpr3()


def ComparisonIdentity_Structural(__unit: None=None) -> IComparer_1[Any]:
    return ComparisonIdentity_FromFunction(compare)


__all__ = ["HashIdentity_FromFunctions", "HashIdentity_Structural", "HashIdentity_Reference", "ComparisonIdentity_FromFunction", "ComparisonIdentity_Structural"]

