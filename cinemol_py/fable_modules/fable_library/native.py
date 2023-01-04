from __future__ import annotations
from abc import abstractmethod
from typing import (TypeVar, Protocol, Generic, Any)
from .types import Array

_T = TypeVar("_T")

class Cons_1(Protocol, Generic[_T]):
    @abstractmethod
    def Allocate(self, len: int) -> Array[_T]:
        ...


def Helpers_allocateArrayFromCons(cons: Cons_1[_T], len_1: int) -> Array[_T]:
    if cons is None:
        return (list)([None]*len_1)

    else: 
        return cons([0]*len_1)



def Helpers_fillImpl(array: Array[_T], value: _T, start: int, count: int) -> Array[_T]:
    for i in range(0, (count - 1) + 1, 1):
        array[i + start] = value
    return array


def Helpers_spliceImpl(array: Array[_T], start: int, delete_count: int) -> Array[_T]:
    for _ in range(1, delete_count + 1, 1):
        array.pop(start)
    return array


def Helpers_indexOfImpl(array: Array[_T], item: _T, start: int) -> Any:
    try: 
        return array.index(item, start)

    except Exception as ex:
        return -1



def Helpers_copyToTypedArray(src: Array[_T], srci: int, trg: Array[_T], trgi: int, cnt: int) -> None:
    diff: int = (trgi - srci) or 0
    for i in range(srci, ((srci + cnt) - 1) + 1, 1):
        trg[i + diff] = src[i]


__all__ = ["Helpers_allocateArrayFromCons", "Helpers_fillImpl", "Helpers_spliceImpl", "Helpers_indexOfImpl", "Helpers_copyToTypedArray"]

