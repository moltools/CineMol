from __future__ import annotations
from dataclasses import dataclass
from typing import (TypeVar, Any, Generic, Optional, Tuple, Callable)
from .array import (fill, map as map_2)
from .list import (FSharpList, cons, empty as empty_1, fold as fold_1, is_empty as is_empty_1, tail, head, of_array_with_tail, singleton)
from .option import (value as value_1, some)
from .reflection import (TypeInfo, class_type, option_type, list_type, bool_type, record_type)
from .seq import (unfold, map as map_1, compare_with, iterate as iterate_1, pick as pick_1, try_pick as try_pick_1)
from .string import (join, format)
from .types import (Array, Record, FSharpRef)
from .util import (IComparer_1, IEnumerator, IEnumerable_1, is_array_like, get_enumerator, equals, to_iterator, compare, ignore, ICollection, structural_hash)

_KEY = TypeVar("_KEY")

_VALUE = TypeVar("_VALUE")

__A = TypeVar("__A")

__B = TypeVar("__B")

__C = TypeVar("__C")

_RESULT = TypeVar("_RESULT")

_A = TypeVar("_A")

_B = TypeVar("_B")

_Z = TypeVar("_Z")

_STATE = TypeVar("_STATE")

_T = TypeVar("_T")

_K = TypeVar("_K")

_V = TypeVar("_V")

def _expr249(gen0: TypeInfo, gen1: TypeInfo) -> TypeInfo:
    return class_type("Map.MapTreeLeaf`2", [gen0, gen1], MapTreeLeaf_2)


class MapTreeLeaf_2(Generic[_KEY, _VALUE]):
    def __init__(self, k: Any, v: Any) -> None:
        self.k: _KEY = k
        self.v: _VALUE = v


MapTreeLeaf_2_reflection = _expr249

def MapTreeLeaf_2__ctor_5BDDA1(k: Any, v: Any) -> MapTreeLeaf_2[_KEY, _VALUE]:
    return MapTreeLeaf_2(k, v)


def MapTreeLeaf_2__get_Key(_: MapTreeLeaf_2[_KEY, Any]) -> _KEY:
    return _.k


def MapTreeLeaf_2__get_Value(_: MapTreeLeaf_2[Any, _VALUE]) -> _VALUE:
    return _.v


def _expr250(gen0: TypeInfo, gen1: TypeInfo) -> TypeInfo:
    return class_type("Map.MapTreeNode`2", [gen0, gen1], MapTreeNode_2, MapTreeLeaf_2_reflection(gen0, gen1))


class MapTreeNode_2(MapTreeLeaf_2, Generic[_KEY, _VALUE]):
    def __init__(self, k: _KEY, v: _VALUE, left: Optional[MapTreeLeaf_2[_KEY, _VALUE]], right: Optional[MapTreeLeaf_2[_KEY, _VALUE]], h: int) -> None:
        super().__init__(k, v)
        self.left: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = left
        self.right: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = right
        self.h: int = h or 0


MapTreeNode_2_reflection = _expr250

def MapTreeNode_2__ctor_Z39DE9543(k: _KEY, v: _VALUE, left: Optional[MapTreeLeaf_2[_KEY, _VALUE]], right: Optional[MapTreeLeaf_2[_KEY, _VALUE]], h: int) -> MapTreeNode_2[_KEY, _VALUE]:
    return MapTreeNode_2(k, v, left, right, h)


def MapTreeNode_2__get_Left(_: MapTreeNode_2[_KEY, _VALUE]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    return _.left


def MapTreeNode_2__get_Right(_: MapTreeNode_2[_KEY, _VALUE]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    return _.right


def MapTreeNode_2__get_Height(_: MapTreeNode_2[Any, Any]) -> int:
    return _.h


def MapTreeModule_empty(__unit: None=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    return None


def MapTreeModule_sizeAux(acc_mut: int, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> int:
    while True:
        (acc, m) = (acc_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                acc_mut = MapTreeModule_sizeAux(acc + 1, MapTreeNode_2__get_Left(m2))
                m_mut = MapTreeNode_2__get_Right(m2)
                continue

            else: 
                return acc + 1


        else: 
            return acc

        break


def MapTreeModule_size(x: Optional[MapTreeLeaf_2[__A, __B]]=None) -> int:
    return MapTreeModule_sizeAux(0, x)


def MapTreeModule_mk(l: Optional[MapTreeLeaf_2[_KEY, _VALUE]], k: _KEY, v: _VALUE, r: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    hl: int
    m: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = l
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        hl = MapTreeNode_2__get_Height(m2) if isinstance(m2, MapTreeNode_2) else 1

    else: 
        hl = 0

    hr: int
    m_1: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = r
    if m_1 is not None:
        m2_1: MapTreeLeaf_2[_KEY, _VALUE] = m_1
        hr = MapTreeNode_2__get_Height(m2_1) if isinstance(m2_1, MapTreeNode_2) else 1

    else: 
        hr = 0

    m_2: int = (hr if (hl < hr) else hl) or 0
    if m_2 == 0:
        return MapTreeLeaf_2__ctor_5BDDA1(k, v)

    else: 
        return MapTreeNode_2__ctor_Z39DE9543(k, v, l, r, m_2 + 1)



def MapTreeModule_rebalance(t1: Optional[MapTreeLeaf_2[_KEY, _VALUE]], k: _KEY, v: _VALUE, t2: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    t1h: int
    m: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = t1
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        t1h = MapTreeNode_2__get_Height(m2) if isinstance(m2, MapTreeNode_2) else 1

    else: 
        t1h = 0

    t2h: int
    m_1: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = t2
    if m_1 is not None:
        m2_1: MapTreeLeaf_2[_KEY, _VALUE] = m_1
        t2h = MapTreeNode_2__get_Height(m2_1) if isinstance(m2_1, MapTreeNode_2) else 1

    else: 
        t2h = 0

    if t2h > (t1h + 2):
        match_value: MapTreeLeaf_2[_KEY, _VALUE] = value_1(t2)
        if isinstance(match_value, MapTreeNode_2):
            def _arrow251(__unit: None=None, t1: Any=t1, k: _KEY=k, v: _VALUE=v, t2: Any=t2) -> int:
                m_2: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = MapTreeNode_2__get_Left(match_value)
                if m_2 is not None:
                    m2_2: MapTreeLeaf_2[_KEY, _VALUE] = m_2
                    return MapTreeNode_2__get_Height(m2_2) if isinstance(m2_2, MapTreeNode_2) else 1

                else: 
                    return 0


            if _arrow251() > (t1h + 1):
                match_value_1: MapTreeLeaf_2[_KEY, _VALUE] = value_1(MapTreeNode_2__get_Left(match_value))
                if isinstance(match_value_1, MapTreeNode_2):
                    return MapTreeModule_mk(MapTreeModule_mk(t1, k, v, MapTreeNode_2__get_Left(match_value_1)), MapTreeLeaf_2__get_Key(match_value_1), MapTreeLeaf_2__get_Value(match_value_1), MapTreeModule_mk(MapTreeNode_2__get_Right(match_value_1), MapTreeLeaf_2__get_Key(match_value), MapTreeLeaf_2__get_Value(match_value), MapTreeNode_2__get_Right(match_value)))

                else: 
                    raise Exception("internal error: Map.rebalance")


            else: 
                return MapTreeModule_mk(MapTreeModule_mk(t1, k, v, MapTreeNode_2__get_Left(match_value)), MapTreeLeaf_2__get_Key(match_value), MapTreeLeaf_2__get_Value(match_value), MapTreeNode_2__get_Right(match_value))


        else: 
            raise Exception("internal error: Map.rebalance")


    elif t1h > (t2h + 2):
        match_value_2: MapTreeLeaf_2[_KEY, _VALUE] = value_1(t1)
        if isinstance(match_value_2, MapTreeNode_2):
            def _arrow252(__unit: None=None, t1: Any=t1, k: _KEY=k, v: _VALUE=v, t2: Any=t2) -> int:
                m_3: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = MapTreeNode_2__get_Right(match_value_2)
                if m_3 is not None:
                    m2_3: MapTreeLeaf_2[_KEY, _VALUE] = m_3
                    return MapTreeNode_2__get_Height(m2_3) if isinstance(m2_3, MapTreeNode_2) else 1

                else: 
                    return 0


            if _arrow252() > (t2h + 1):
                match_value_3: MapTreeLeaf_2[_KEY, _VALUE] = value_1(MapTreeNode_2__get_Right(match_value_2))
                if isinstance(match_value_3, MapTreeNode_2):
                    return MapTreeModule_mk(MapTreeModule_mk(MapTreeNode_2__get_Left(match_value_2), MapTreeLeaf_2__get_Key(match_value_2), MapTreeLeaf_2__get_Value(match_value_2), MapTreeNode_2__get_Left(match_value_3)), MapTreeLeaf_2__get_Key(match_value_3), MapTreeLeaf_2__get_Value(match_value_3), MapTreeModule_mk(MapTreeNode_2__get_Right(match_value_3), k, v, t2))

                else: 
                    raise Exception("internal error: Map.rebalance")


            else: 
                return MapTreeModule_mk(MapTreeNode_2__get_Left(match_value_2), MapTreeLeaf_2__get_Key(match_value_2), MapTreeLeaf_2__get_Value(match_value_2), MapTreeModule_mk(MapTreeNode_2__get_Right(match_value_2), k, v, t2))


        else: 
            raise Exception("internal error: Map.rebalance")


    else: 
        return MapTreeModule_mk(t1, k, v, t2)



def MapTreeModule_add(comparer: IComparer_1[_KEY], k: _KEY, v: _VALUE, m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        c: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
        if isinstance(m2, MapTreeNode_2):
            if c < 0:
                return MapTreeModule_rebalance(MapTreeModule_add(comparer, k, v, MapTreeNode_2__get_Left(m2)), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeNode_2__get_Right(m2))

            elif c == 0:
                return MapTreeNode_2__ctor_Z39DE9543(k, v, MapTreeNode_2__get_Left(m2), MapTreeNode_2__get_Right(m2), MapTreeNode_2__get_Height(m2))

            else: 
                return MapTreeModule_rebalance(MapTreeNode_2__get_Left(m2), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeModule_add(comparer, k, v, MapTreeNode_2__get_Right(m2)))


        elif c < 0:
            return MapTreeNode_2__ctor_Z39DE9543(k, v, MapTreeModule_empty(), m, 2)

        elif c == 0:
            return MapTreeLeaf_2__ctor_5BDDA1(k, v)

        else: 
            return MapTreeNode_2__ctor_Z39DE9543(k, v, m, MapTreeModule_empty(), 2)


    else: 
        return MapTreeLeaf_2__ctor_5BDDA1(k, v)



def MapTreeModule_tryFind(comparer_mut: IComparer_1[_KEY], k_mut: _KEY, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> Optional[_VALUE]:
    while True:
        (comparer, k, m) = (comparer_mut, k_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            c: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
            if c == 0:
                return some(MapTreeLeaf_2__get_Value(m2))

            elif isinstance(m2, MapTreeNode_2):
                comparer_mut = comparer
                k_mut = k
                m_mut = MapTreeNode_2__get_Left(m2) if (c < 0) else MapTreeNode_2__get_Right(m2)
                continue

            else: 
                return None


        else: 
            return None

        break


def MapTreeModule_find(comparer: IComparer_1[_KEY], k: _KEY, m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> _VALUE:
    match_value: Optional[_VALUE] = MapTreeModule_tryFind(comparer, k, m)
    if match_value is None:
        raise Exception()

    else: 
        return value_1(match_value)



def MapTreeModule_partition1(comparer: IComparer_1[_KEY], f: Any, k: _KEY, v: __A, acc1: Optional[MapTreeLeaf_2[_KEY, __A]]=None, acc2: Optional[MapTreeLeaf_2[_KEY, __A]]=None) -> Tuple[Optional[MapTreeLeaf_2[_KEY, __A]], Optional[MapTreeLeaf_2[_KEY, __A]]]:
    if f(k, v):
        return (MapTreeModule_add(comparer, k, v, acc1), acc2)

    else: 
        return (acc1, MapTreeModule_add(comparer, k, v, acc2))



def MapTreeModule_partitionAux(comparer_mut: IComparer_1[_KEY], f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], acc__mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], acc__1_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> Tuple[Optional[MapTreeLeaf_2[_KEY, _VALUE]], Optional[MapTreeLeaf_2[_KEY, _VALUE]]]:
    while True:
        (comparer, f, m, acc_, acc__1) = (comparer_mut, f_mut, m_mut, acc__mut, acc__1_mut)
        acc: Tuple[Optional[MapTreeLeaf_2[_KEY, _VALUE]], Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = (acc_, acc__1)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                acc_1: Tuple[Optional[MapTreeLeaf_2[_KEY, _VALUE]], Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_partitionAux(comparer, f, MapTreeNode_2__get_Right(m2), acc[0], acc[1])
                acc_4: Tuple[Optional[MapTreeLeaf_2[_KEY, _VALUE]], Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_partition1(comparer, f, MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), acc_1[0], acc_1[1])
                comparer_mut = comparer
                f_mut = f
                m_mut = MapTreeNode_2__get_Left(m2)
                acc__mut = acc_4[0]
                acc__1_mut = acc_4[1]
                continue

            else: 
                return MapTreeModule_partition1(comparer, f, MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), acc[0], acc[1])


        else: 
            return acc

        break


def MapTreeModule_partition(comparer: IComparer_1[_KEY], f: Callable[[_KEY, __A], bool], m: Optional[MapTreeLeaf_2[_KEY, __A]]=None) -> Tuple[Optional[MapTreeLeaf_2[_KEY, __A]], Optional[MapTreeLeaf_2[_KEY, __A]]]:
    return MapTreeModule_partitionAux(comparer, f, m, MapTreeModule_empty(), MapTreeModule_empty())


def MapTreeModule_filter1(comparer: IComparer_1[_KEY], f: Any, k: _KEY, v: __A, acc: Optional[MapTreeLeaf_2[_KEY, __A]]=None) -> Optional[MapTreeLeaf_2[_KEY, __A]]:
    if f(k, v):
        return MapTreeModule_add(comparer, k, v, acc)

    else: 
        return acc



def MapTreeModule_filterAux(comparer_mut: IComparer_1[_KEY], f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], acc_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    while True:
        (comparer, f, m, acc) = (comparer_mut, f_mut, m_mut, acc_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                acc_1: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = MapTreeModule_filterAux(comparer, f, MapTreeNode_2__get_Left(m2), acc)
                acc_2: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = MapTreeModule_filter1(comparer, f, MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), acc_1)
                comparer_mut = comparer
                f_mut = f
                m_mut = MapTreeNode_2__get_Right(m2)
                acc_mut = acc_2
                continue

            else: 
                return MapTreeModule_filter1(comparer, f, MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), acc)


        else: 
            return acc

        break


def MapTreeModule_filter(comparer: IComparer_1[_KEY], f: Callable[[_KEY, __A], bool], m: Optional[MapTreeLeaf_2[_KEY, __A]]=None) -> Optional[MapTreeLeaf_2[_KEY, __A]]:
    return MapTreeModule_filterAux(comparer, f, m, MapTreeModule_empty())


def MapTreeModule_spliceOutSuccessor(m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Tuple[_KEY, _VALUE, Optional[MapTreeLeaf_2[_KEY, _VALUE]]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        if isinstance(m2, MapTreeNode_2):
            if MapTreeNode_2__get_Left(m2) is None:
                return (MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeNode_2__get_Right(m2))

            else: 
                pattern_input: Tuple[_KEY, _VALUE, Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_spliceOutSuccessor(MapTreeNode_2__get_Left(m2))
                return (pattern_input[0], pattern_input[1], MapTreeModule_mk(pattern_input[2], MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeNode_2__get_Right(m2)))


        else: 
            return (MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeModule_empty())


    else: 
        raise Exception("internal error: Map.spliceOutSuccessor")



def MapTreeModule_remove(comparer: IComparer_1[_KEY], k: _KEY, m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        c: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
        if isinstance(m2, MapTreeNode_2):
            if c < 0:
                return MapTreeModule_rebalance(MapTreeModule_remove(comparer, k, MapTreeNode_2__get_Left(m2)), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeNode_2__get_Right(m2))

            elif c == 0:
                if MapTreeNode_2__get_Left(m2) is None:
                    return MapTreeNode_2__get_Right(m2)

                elif MapTreeNode_2__get_Right(m2) is None:
                    return MapTreeNode_2__get_Left(m2)

                else: 
                    pattern_input: Tuple[_KEY, _VALUE, Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_spliceOutSuccessor(MapTreeNode_2__get_Right(m2))
                    return MapTreeModule_mk(MapTreeNode_2__get_Left(m2), pattern_input[0], pattern_input[1], pattern_input[2])


            else: 
                return MapTreeModule_rebalance(MapTreeNode_2__get_Left(m2), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeModule_remove(comparer, k, MapTreeNode_2__get_Right(m2)))


        elif c == 0:
            return MapTreeModule_empty()

        else: 
            return m


    else: 
        return MapTreeModule_empty()



def MapTreeModule_change(comparer: IComparer_1[_KEY], k: _KEY, u: Callable[[Optional[_VALUE]], Optional[_VALUE]], m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        if isinstance(m2, MapTreeNode_2):
            c: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
            if c < 0:
                return MapTreeModule_rebalance(MapTreeModule_change(comparer, k, u, MapTreeNode_2__get_Left(m2)), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeNode_2__get_Right(m2))

            elif c == 0:
                match_value_1: Optional[_VALUE] = u(some(MapTreeLeaf_2__get_Value(m2)))
                if match_value_1 is not None:
                    return MapTreeNode_2__ctor_Z39DE9543(k, value_1(match_value_1), MapTreeNode_2__get_Left(m2), MapTreeNode_2__get_Right(m2), MapTreeNode_2__get_Height(m2))

                elif MapTreeNode_2__get_Left(m2) is None:
                    return MapTreeNode_2__get_Right(m2)

                elif MapTreeNode_2__get_Right(m2) is None:
                    return MapTreeNode_2__get_Left(m2)

                else: 
                    pattern_input: Tuple[_KEY, _VALUE, Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_spliceOutSuccessor(MapTreeNode_2__get_Right(m2))
                    return MapTreeModule_mk(MapTreeNode_2__get_Left(m2), pattern_input[0], pattern_input[1], pattern_input[2])


            else: 
                return MapTreeModule_rebalance(MapTreeNode_2__get_Left(m2), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), MapTreeModule_change(comparer, k, u, MapTreeNode_2__get_Right(m2)))


        else: 
            c_1: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
            if c_1 < 0:
                match_value_2: Optional[_VALUE] = u(None)
                if match_value_2 is not None:
                    return MapTreeNode_2__ctor_Z39DE9543(k, value_1(match_value_2), MapTreeModule_empty(), m, 2)

                else: 
                    return m


            elif c_1 == 0:
                match_value_3: Optional[_VALUE] = u(some(MapTreeLeaf_2__get_Value(m2)))
                if match_value_3 is not None:
                    return MapTreeLeaf_2__ctor_5BDDA1(k, value_1(match_value_3))

                else: 
                    return MapTreeModule_empty()


            else: 
                match_value_4: Optional[_VALUE] = u(None)
                if match_value_4 is not None:
                    return MapTreeNode_2__ctor_Z39DE9543(k, value_1(match_value_4), m, MapTreeModule_empty(), 2)

                else: 
                    return m




    else: 
        match_value: Optional[_VALUE] = u(None)
        if match_value is not None:
            return MapTreeLeaf_2__ctor_5BDDA1(k, value_1(match_value))

        else: 
            return m




def MapTreeModule_mem(comparer_mut: IComparer_1[_KEY], k_mut: _KEY, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> bool:
    while True:
        (comparer, k, m) = (comparer_mut, k_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            c: int = comparer.Compare(k, MapTreeLeaf_2__get_Key(m2)) or 0
            if isinstance(m2, MapTreeNode_2):
                if c < 0:
                    comparer_mut = comparer
                    k_mut = k
                    m_mut = MapTreeNode_2__get_Left(m2)
                    continue

                elif c == 0:
                    return True

                else: 
                    comparer_mut = comparer
                    k_mut = k
                    m_mut = MapTreeNode_2__get_Right(m2)
                    continue


            else: 
                return c == 0


        else: 
            return False

        break


def MapTreeModule_iterOpt(f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> None:
    while True:
        (f, m) = (f_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                MapTreeModule_iterOpt(f, MapTreeNode_2__get_Left(m2))
                f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))
                f_mut = f
                m_mut = MapTreeNode_2__get_Right(m2)
                continue

            else: 
                f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))


        break


def MapTreeModule_iter(f: Callable[[__A, __B], None], m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> None:
    MapTreeModule_iterOpt(f, m)


def MapTreeModule_tryPickOpt(f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> Optional[__A]:
    while True:
        (f, m) = (f_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                match_value: Optional[__A] = MapTreeModule_tryPickOpt(f, MapTreeNode_2__get_Left(m2))
                if match_value is None:
                    match_value_1: Optional[__A] = f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))
                    if match_value_1 is None:
                        f_mut = f
                        m_mut = MapTreeNode_2__get_Right(m2)
                        continue

                    else: 
                        return match_value_1


                else: 
                    return match_value


            else: 
                return f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))


        else: 
            return None

        break


def MapTreeModule_tryPick(f: Callable[[__A, __B], Optional[__C]], m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> Optional[__C]:
    return MapTreeModule_tryPickOpt(f, m)


def MapTreeModule_existsOpt(f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> bool:
    while True:
        (f, m) = (f_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                if True if MapTreeModule_existsOpt(f, MapTreeNode_2__get_Left(m2)) else f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)):
                    return True

                else: 
                    f_mut = f
                    m_mut = MapTreeNode_2__get_Right(m2)
                    continue


            else: 
                return f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))


        else: 
            return False

        break


def MapTreeModule_exists(f: Callable[[__A, __B], bool], m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> bool:
    return MapTreeModule_existsOpt(f, m)


def MapTreeModule_forallOpt(f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> bool:
    while True:
        (f, m) = (f_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                if f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)) if MapTreeModule_forallOpt(f, MapTreeNode_2__get_Left(m2)) else False:
                    f_mut = f
                    m_mut = MapTreeNode_2__get_Right(m2)
                    continue

                else: 
                    return False


            else: 
                return f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))


        else: 
            return True

        break


def MapTreeModule_forall(f: Callable[[__A, __B], bool], m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> bool:
    return MapTreeModule_forallOpt(f, m)


def MapTreeModule_map(f: Callable[[_VALUE], _RESULT], m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _RESULT]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        if isinstance(m2, MapTreeNode_2):
            l2: Optional[MapTreeLeaf_2[_KEY, _RESULT]] = MapTreeModule_map(f, MapTreeNode_2__get_Left(m2))
            v2: _RESULT = f(MapTreeLeaf_2__get_Value(m2))
            r2: Optional[MapTreeLeaf_2[_KEY, _RESULT]] = MapTreeModule_map(f, MapTreeNode_2__get_Right(m2))
            return MapTreeNode_2__ctor_Z39DE9543(MapTreeLeaf_2__get_Key(m2), v2, l2, r2, MapTreeNode_2__get_Height(m2))

        else: 
            return MapTreeLeaf_2__ctor_5BDDA1(MapTreeLeaf_2__get_Key(m2), f(MapTreeLeaf_2__get_Value(m2)))


    else: 
        return MapTreeModule_empty()



def MapTreeModule_mapiOpt(f: Any, m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> Optional[MapTreeLeaf_2[_KEY, _RESULT]]:
    if m is not None:
        m2: MapTreeLeaf_2[_KEY, _VALUE] = m
        if isinstance(m2, MapTreeNode_2):
            l2: Optional[MapTreeLeaf_2[_KEY, _RESULT]] = MapTreeModule_mapiOpt(f, MapTreeNode_2__get_Left(m2))
            v2: _RESULT = f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))
            r2: Optional[MapTreeLeaf_2[_KEY, _RESULT]] = MapTreeModule_mapiOpt(f, MapTreeNode_2__get_Right(m2))
            return MapTreeNode_2__ctor_Z39DE9543(MapTreeLeaf_2__get_Key(m2), v2, l2, r2, MapTreeNode_2__get_Height(m2))

        else: 
            return MapTreeLeaf_2__ctor_5BDDA1(MapTreeLeaf_2__get_Key(m2), f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)))


    else: 
        return MapTreeModule_empty()



def MapTreeModule_mapi(f: Callable[[__A, __B], __C], m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> Optional[MapTreeLeaf_2[__A, __C]]:
    return MapTreeModule_mapiOpt(f, m)


def MapTreeModule_foldBackOpt(f_mut: Any, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], x_mut: __A) -> __A:
    while True:
        (f, m, x) = (f_mut, m_mut, x_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                x_1: __A = MapTreeModule_foldBackOpt(f, MapTreeNode_2__get_Right(m2), x)
                x_2: __A = f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), x_1)
                f_mut = f
                m_mut = MapTreeNode_2__get_Left(m2)
                x_mut = x_2
                continue

            else: 
                return f(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), x)


        else: 
            return x

        break


def MapTreeModule_foldBack(f: Callable[[__A, __B, __C], __C], m: Optional[MapTreeLeaf_2[__A, __B]], x: __C) -> __C:
    return MapTreeModule_foldBackOpt(f, m, x)


def MapTreeModule_foldOpt(f_mut: Any, x_mut: __A, m_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]]) -> __A:
    while True:
        (f, x, m) = (f_mut, x_mut, m_mut)
        if m is not None:
            m2: MapTreeLeaf_2[_KEY, _VALUE] = m
            if isinstance(m2, MapTreeNode_2):
                f_mut = f
                x_mut = f(MapTreeModule_foldOpt(f, x, MapTreeNode_2__get_Left(m2)), MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))
                m_mut = MapTreeNode_2__get_Right(m2)
                continue

            else: 
                return f(x, MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2))


        else: 
            return x

        break


def MapTreeModule_fold(f: Callable[[__A, __B, __C], __A], x: __A, m: Optional[MapTreeLeaf_2[__B, __C]]=None) -> __A:
    return MapTreeModule_foldOpt(f, x, m)


def MapTreeModule_foldSectionOpt(comparer: IComparer_1[_KEY], lo: _KEY, hi: _KEY, f: Any, m: Optional[MapTreeLeaf_2[_KEY, _VALUE]], x: _A) -> _A:
    def fold_from_to(f_1_mut: Any, m_1_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], x_1_mut: _A, comparer: IComparer_1[Any]=comparer, lo: _KEY=lo, hi: _KEY=hi, f: Any=f, m: Any=m, x: _A=x) -> _A:
        while True:
            (f_1, m_1, x_1) = (f_1_mut, m_1_mut, x_1_mut)
            if m_1 is not None:
                m2: MapTreeLeaf_2[_KEY, _VALUE] = m_1
                if isinstance(m2, MapTreeNode_2):
                    c_lo_key: int = comparer.Compare(lo, MapTreeLeaf_2__get_Key(m2)) or 0
                    c_key_hi: int = comparer.Compare(MapTreeLeaf_2__get_Key(m2), hi) or 0
                    x_2: _A = fold_from_to(f_1, MapTreeNode_2__get_Left(m2), x_1) if (c_lo_key < 0) else x_1
                    x_3: _A = f_1(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), x_2) if ((c_key_hi <= 0) if (c_lo_key <= 0) else False) else x_2
                    if c_key_hi < 0:
                        f_1_mut = f_1
                        m_1_mut = MapTreeNode_2__get_Right(m2)
                        x_1_mut = x_3
                        continue

                    else: 
                        return x_3


                elif (comparer.Compare(MapTreeLeaf_2__get_Key(m2), hi) <= 0) if (comparer.Compare(lo, MapTreeLeaf_2__get_Key(m2)) <= 0) else False:
                    return f_1(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2), x_1)

                else: 
                    return x_1


            else: 
                return x_1

            break

    if comparer.Compare(lo, hi) == 1:
        return x

    else: 
        return fold_from_to(f, m, x)



def MapTreeModule_foldSection(comparer: IComparer_1[_KEY], lo: _KEY, hi: _KEY, f: Callable[[_KEY, __A, __B], __B], m: Optional[MapTreeLeaf_2[_KEY, __A]], x: __B) -> __B:
    return MapTreeModule_foldSectionOpt(comparer, lo, hi, f, m, x)


def MapTreeModule_toList(m: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> FSharpList[Tuple[_KEY, _VALUE]]:
    def loop(m_1_mut: Optional[MapTreeLeaf_2[_KEY, _VALUE]], acc_mut: FSharpList[Tuple[_KEY, _VALUE]], m: Any=m) -> FSharpList[Tuple[_KEY, _VALUE]]:
        while True:
            (m_1, acc) = (m_1_mut, acc_mut)
            if m_1 is not None:
                m2: MapTreeLeaf_2[_KEY, _VALUE] = m_1
                if isinstance(m2, MapTreeNode_2):
                    m_1_mut = MapTreeNode_2__get_Left(m2)
                    acc_mut = cons((MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)), loop(MapTreeNode_2__get_Right(m2), acc))
                    continue

                else: 
                    return cons((MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)), acc)


            else: 
                return acc

            break

    return loop(m, empty_1())


def MapTreeModule_copyToArray(m: Optional[MapTreeLeaf_2[__A, __B]], arr: Array[Any], i: int) -> None:
    j: int = i or 0
    def _arrow253(x: __A, y: __B, m: Any=m, arr: Any=arr, i: int=i) -> None:
        nonlocal j
        arr[j] = (x, y)
        j = (j + 1) or 0

    MapTreeModule_iter(_arrow253, m)


def MapTreeModule_toArray(m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> Array[Any]:
    n: int = MapTreeModule_size(m) or 0
    res: Array[Any] = fill([0] * n, 0, n, (None, None))
    MapTreeModule_copyToArray(m, res, 0)
    return res


def MapTreeModule_ofList(comparer: IComparer_1[__A], l: FSharpList[Tuple[__A, __B]]) -> Optional[MapTreeLeaf_2[__A, __B]]:
    def _arrow254(acc: Optional[MapTreeLeaf_2[__A, __B]], tupled_arg: Tuple[__A, __B], comparer: IComparer_1[Any]=comparer, l: Any=l) -> Optional[MapTreeLeaf_2[__A, __B]]:
        return MapTreeModule_add(comparer, tupled_arg[0], tupled_arg[1], acc)

    return fold_1(_arrow254, MapTreeModule_empty(), l)


def MapTreeModule_mkFromEnumerator(comparer_mut: IComparer_1[__A], acc_mut: Optional[MapTreeLeaf_2[__A, __B]], e_mut: IEnumerator[Tuple[__A, __B]]) -> Optional[MapTreeLeaf_2[__A, __B]]:
    while True:
        (comparer, acc, e) = (comparer_mut, acc_mut, e_mut)
        if e.System_Collections_IEnumerator_MoveNext():
            pattern_input: Tuple[__A, __B] = e.System_Collections_Generic_IEnumerator_1_get_Current()
            comparer_mut = comparer
            acc_mut = MapTreeModule_add(comparer, pattern_input[0], pattern_input[1], acc)
            e_mut = e
            continue

        else: 
            return acc

        break


def MapTreeModule_ofArray(comparer: IComparer_1[_KEY], arr: Array[Tuple[_KEY, _VALUE]]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    res: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = MapTreeModule_empty()
    for idx in range(0, (len(arr) - 1) + 1, 1):
        for_loop_var: Tuple[_KEY, _VALUE] = arr[idx]
        res = MapTreeModule_add(comparer, for_loop_var[0], for_loop_var[1], res)
    return res


def MapTreeModule_ofSeq(comparer: IComparer_1[_KEY], c: IEnumerable_1[Tuple[_KEY, _VALUE]]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    if is_array_like(c):
        return MapTreeModule_ofArray(comparer, c)

    elif isinstance(c, FSharpList):
        return MapTreeModule_ofList(comparer, c)

    else: 
        with get_enumerator(c) as ie:
            return MapTreeModule_mkFromEnumerator(comparer, MapTreeModule_empty(), ie)



def _expr255(gen0: TypeInfo, gen1: TypeInfo) -> TypeInfo:
    return record_type("Map.MapTreeModule.MapIterator`2", [gen0, gen1], MapTreeModule_MapIterator_2, lambda: [("stack", list_type(option_type(MapTreeLeaf_2_reflection(gen0, gen1)))), ("started", bool_type)])


@dataclass(eq = False, repr = False)
class MapTreeModule_MapIterator_2(Record, Generic[_KEY, _VALUE]):
    stack: FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]]
    started: bool

MapTreeModule_MapIterator_2_reflection = _expr255

def MapTreeModule_collapseLHS(stack_mut: FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]]) -> FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]]:
    while True:
        (stack,) = (stack_mut,)
        if not is_empty_1(stack):
            rest: FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = tail(stack)
            m: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = head(stack)
            if m is not None:
                m2: MapTreeLeaf_2[_KEY, _VALUE] = m
                if isinstance(m2, MapTreeNode_2):
                    stack_mut = of_array_with_tail([MapTreeNode_2__get_Left(m2), MapTreeLeaf_2__ctor_5BDDA1(MapTreeLeaf_2__get_Key(m2), MapTreeLeaf_2__get_Value(m2)), MapTreeNode_2__get_Right(m2)], rest)
                    continue

                else: 
                    return stack


            else: 
                stack_mut = rest
                continue


        else: 
            return empty_1()

        break


def MapTreeModule_mkIterator(m: Optional[MapTreeLeaf_2[__A, __B]]=None) -> MapTreeModule_MapIterator_2[__A, __B]:
    return MapTreeModule_MapIterator_2(MapTreeModule_collapseLHS(singleton(m)), False)


def MapTreeModule_notStarted(__unit: None=None) -> Any:
    raise Exception("enumeration not started")


def MapTreeModule_alreadyFinished(__unit: None=None) -> Any:
    raise Exception("enumeration already finished")


def MapTreeModule_current(i: MapTreeModule_MapIterator_2[_KEY, _VALUE]) -> Any:
    if i.started:
        match_value: FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = i.stack
        if not is_empty_1(match_value):
            if head(match_value) is not None:
                m: MapTreeLeaf_2[_KEY, _VALUE] = head(match_value)
                if isinstance(m, MapTreeNode_2):
                    raise Exception("Please report error: Map iterator, unexpected stack for current")

                else: 
                    return (MapTreeLeaf_2__get_Key(m), MapTreeLeaf_2__get_Value(m))


            else: 
                raise Exception("Please report error: Map iterator, unexpected stack for current")


        else: 
            return MapTreeModule_alreadyFinished()


    else: 
        return MapTreeModule_notStarted()



def MapTreeModule_moveNext(i: MapTreeModule_MapIterator_2[Any, Any]) -> bool:
    if i.started:
        match_value: FSharpList[Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = i.stack
        if not is_empty_1(match_value):
            if head(match_value) is not None:
                m: MapTreeLeaf_2[_KEY, _VALUE] = head(match_value)
                if isinstance(m, MapTreeNode_2):
                    raise Exception("Please report error: Map iterator, unexpected stack for moveNext")

                else: 
                    i.stack = MapTreeModule_collapseLHS(tail(match_value))
                    return not is_empty_1(i.stack)


            else: 
                raise Exception("Please report error: Map iterator, unexpected stack for moveNext")


        else: 
            return False


    else: 
        i.started = True
        return not is_empty_1(i.stack)



def MapTreeModule_mkIEnumerator(m: Optional[MapTreeLeaf_2[_A, _B]]=None) -> IEnumerator[Any]:
    i: MapTreeModule_MapIterator_2[_A, _B] = MapTreeModule_mkIterator(m)
    class ObjectExpr256(IEnumerator[Any]):
        def System_Collections_Generic_IEnumerator_1_get_Current(self, __unit: None=None, m: Any=m) -> Any:
            return MapTreeModule_current(i)

        def System_Collections_IEnumerator_get_Current(self, __unit: None=None, m: Any=m) -> Any:
            return MapTreeModule_current(i)

        def System_Collections_IEnumerator_MoveNext(self, __unit: None=None, m: Any=m) -> bool:
            return MapTreeModule_moveNext(i)

        def System_Collections_IEnumerator_Reset(self, __unit: None=None, m: Any=m) -> None:
            nonlocal i
            i = MapTreeModule_mkIterator(m)

        def Dispose(self, __unit: None=None, m: Any=m) -> None:
            pass

    return ObjectExpr256()


def MapTreeModule_toSeq(s: Optional[MapTreeLeaf_2[__A, __B]]=None) -> IEnumerable_1[Any]:
    def generator(en_1: IEnumerator[Any], s: Any=s) -> Optional[Tuple[Any, IEnumerator[Any]]]:
        if en_1.System_Collections_IEnumerator_MoveNext():
            return (en_1.System_Collections_Generic_IEnumerator_1_get_Current(), en_1)

        else: 
            return None


    return unfold(generator, MapTreeModule_mkIEnumerator(s))


def _expr259(gen0: TypeInfo, gen1: TypeInfo) -> TypeInfo:
    return class_type("Map.FSharpMap", [gen0, gen1], FSharpMap)


class FSharpMap(Generic[_KEY, _VALUE]):
    def __init__(self, comparer: IComparer_1[_KEY], tree: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> None:
        self.comparer: IComparer_1[_KEY] = comparer
        self.tree: Optional[MapTreeLeaf_2[_KEY, _VALUE]] = tree

    def GetHashCode(self, __unit: None=None) -> int:
        this: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__ComputeHashCode(this)

    def __eq__(self, that: Any=None) -> bool:
        this: FSharpMap[_KEY, _VALUE] = self
        if isinstance(that, FSharpMap):
            with get_enumerator(this) as e1:
                with get_enumerator(that) as e2:
                    def loop(__unit: None=None) -> bool:
                        m1: bool = e1.System_Collections_IEnumerator_MoveNext()
                        if m1 == e2.System_Collections_IEnumerator_MoveNext():
                            if not m1:
                                return True

                            else: 
                                e1c: Any = e1.System_Collections_Generic_IEnumerator_1_get_Current()
                                e2c: Any = e2.System_Collections_Generic_IEnumerator_1_get_Current()
                                if equals(e1c[1], e2c[1]) if equals(e1c[0], e2c[0]) else False:
                                    return loop()

                                else: 
                                    return False



                        else: 
                            return False


                    return loop()

        else: 
            return False


    def __str__(self, __unit: None=None) -> str:
        this: FSharpMap[_KEY, _VALUE] = self
        def _arrow257(kv: Any) -> str:
            return format("({0}, {1})", kv[0], kv[1])

        return ("map [" + join("; ", map_1(_arrow257, this))) + "]"

    @property
    def Symbol_toStringTag(self, __unit: None=None) -> str:
        return "FSharpMap"

    def to_json(self, __unit: None=None) -> Any:
        this: FSharpMap[_KEY, _VALUE] = self
        return Array.from_(this)

    def GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        _: FSharpMap[_KEY, _VALUE] = self
        return MapTreeModule_mkIEnumerator(_.tree)

    def __iter__(self) -> IEnumerator[Any]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        _: FSharpMap[_KEY, _VALUE] = self
        return MapTreeModule_mkIEnumerator(_.tree)

    def __cmp__(self, obj: Any=None) -> int:
        m: FSharpMap[_KEY, _VALUE] = self
        if isinstance(obj, FSharpMap):
            def _arrow258(kvp1: Any, kvp2: Any) -> int:
                c: int = m.comparer.Compare(kvp1[0], kvp2[0]) or 0
                return c if (c != 0) else compare(kvp1[1], kvp2[1])

            return compare_with(_arrow258, m, obj)

        else: 
            raise Exception("not comparable\\nParameter name: obj")


    def System_Collections_Generic_ICollection_1_Add2B595(self, x: Any) -> None:
        ignore(x)
        raise Exception("Map cannot be mutated")

    def System_Collections_Generic_ICollection_1_Clear(self, __unit: None=None) -> None:
        raise Exception("Map cannot be mutated")

    def System_Collections_Generic_ICollection_1_Remove2B595(self, x: Any) -> bool:
        ignore(x)
        raise Exception("Map cannot be mutated")

    def System_Collections_Generic_ICollection_1_Contains2B595(self, x: Any) -> bool:
        m: FSharpMap[_KEY, _VALUE] = self
        return equals(FSharpMap__get_Item(m, x[0]), x[1]) if FSharpMap__ContainsKey(m, x[0]) else False

    def System_Collections_Generic_ICollection_1_CopyToZ3B4C077E(self, arr: Array[Any], i: int) -> None:
        m: FSharpMap[_KEY, _VALUE] = self
        MapTreeModule_copyToArray(m.tree, arr, i)

    def System_Collections_Generic_ICollection_1_get_IsReadOnly(self, __unit: None=None) -> bool:
        return True

    def __len__(self, __unit: None=None) -> int:
        m: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__get_Count(m)

    def __len__(self, __unit: None=None) -> int:
        m: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__get_Count(m)

    @property
    def size(self, __unit: None=None) -> int:
        m: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__get_Count(m)

    def clear(self, __unit: None=None) -> None:
        raise Exception("Map cannot be mutated")

    def delete(self, _arg: Optional[_KEY]=None) -> bool:
        raise Exception("Map cannot be mutated")
        return False

    def entries(self, __unit: None=None) -> IEnumerable_1[Tuple[_KEY, _VALUE]]:
        m: FSharpMap[_KEY, _VALUE] = self
        def mapping(p: Any) -> Tuple[_KEY, _VALUE]:
            return (p[0], p[1])

        return map_1(mapping, m)

    def __getitem__(self, k: Optional[_KEY]=None) -> _VALUE:
        m: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__get_Item(m, k)

    def __contains__(self, k: Optional[_KEY]=None) -> bool:
        m: FSharpMap[_KEY, _VALUE] = self
        return FSharpMap__ContainsKey(m, k)

    def keys(self, __unit: None=None) -> IEnumerable_1[_KEY]:
        m: FSharpMap[_KEY, _VALUE] = self
        def mapping(p: Any) -> _KEY:
            return p[0]

        return map_1(mapping, m)

    def __setitem__(self, k: _KEY, v: _VALUE) -> Map_2[_KEY, _VALUE]:
        m: FSharpMap[_KEY, _VALUE] = self
        raise Exception("Map cannot be mutated")
        return m

    def values(self, __unit: None=None) -> IEnumerable_1[_VALUE]:
        m: FSharpMap[_KEY, _VALUE] = self
        def mapping(p: Any) -> _VALUE:
            return p[1]

        return map_1(mapping, m)

    def for_each(self, f: Callable[[_VALUE, _KEY, Map_2[_KEY, _VALUE]], None], this_arg: Optional[Any]=None) -> None:
        m: FSharpMap[_KEY, _VALUE] = self
        def action(p: Any) -> None:
            f(p[1], p[0], m)

        iterate_1(action, m)


FSharpMap_reflection = _expr259

def FSharpMap__ctor(comparer: IComparer_1[_KEY], tree: Optional[MapTreeLeaf_2[_KEY, _VALUE]]=None) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap(comparer, tree)


def FSharpMap_Empty(comparer: IComparer_1[_KEY]) -> FSharpMap[_KEY, Any]:
    return FSharpMap__ctor(comparer, MapTreeModule_empty())


def FSharpMap__get_Comparer(m: FSharpMap[_KEY, Any]) -> IComparer_1[_KEY]:
    return m.comparer


def FSharpMap__get_Tree(m: FSharpMap[_KEY, _VALUE]) -> Optional[MapTreeLeaf_2[_KEY, _VALUE]]:
    return m.tree


def FSharpMap__Add(m: FSharpMap[_KEY, _VALUE], key: _KEY, value: _VALUE) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_add(m.comparer, key, value, m.tree))


def FSharpMap__Change(m: FSharpMap[_KEY, _VALUE], key: _KEY, f: Callable[[Optional[_VALUE]], Optional[_VALUE]]) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_change(m.comparer, key, f, m.tree))


def FSharpMap__get_IsEmpty(m: FSharpMap[Any, Any]) -> bool:
    return m.tree is None


def FSharpMap__get_Item(m: FSharpMap[_KEY, _VALUE], key: _KEY) -> _VALUE:
    return MapTreeModule_find(m.comparer, key, m.tree)


def FSharpMap__TryPick(m: FSharpMap[_KEY, _VALUE], f: Callable[[_KEY, _VALUE], Optional[__A]]) -> Optional[__A]:
    return MapTreeModule_tryPick(f, m.tree)


def FSharpMap__Exists(m: FSharpMap[_KEY, _VALUE], predicate: Callable[[_KEY, _VALUE], bool]) -> bool:
    return MapTreeModule_exists(predicate, m.tree)


def FSharpMap__Filter(m: FSharpMap[_KEY, _VALUE], predicate: Callable[[_KEY, _VALUE], bool]) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_filter(m.comparer, predicate, m.tree))


def FSharpMap__ForAll(m: FSharpMap[_KEY, _VALUE], predicate: Callable[[_KEY, _VALUE], bool]) -> bool:
    return MapTreeModule_forall(predicate, m.tree)


def FSharpMap__Fold(m: FSharpMap[_KEY, _VALUE], f: Callable[[_KEY, _VALUE, __A], __A], acc: __A) -> __A:
    return MapTreeModule_foldBack(f, m.tree, acc)


def FSharpMap__FoldSection(m: FSharpMap[_KEY, _VALUE], lo: _KEY, hi: _KEY, f: Callable[[_KEY, _VALUE, _Z], _Z], acc: _Z) -> _Z:
    return MapTreeModule_foldSection(m.comparer, lo, hi, f, m.tree, acc)


def FSharpMap__Iterate(m: FSharpMap[_KEY, _VALUE], f: Callable[[_KEY, _VALUE], None]) -> None:
    MapTreeModule_iter(f, m.tree)


def FSharpMap__MapRange(m: FSharpMap[_KEY, _VALUE], f: Callable[[_VALUE], _RESULT]) -> FSharpMap[_KEY, _RESULT]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_map(f, m.tree))


def FSharpMap__Map(m: FSharpMap[_KEY, _VALUE], f: Callable[[_KEY, _VALUE], _B]) -> FSharpMap[_KEY, _B]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_mapi(f, m.tree))


def FSharpMap__Partition(m: FSharpMap[_KEY, _VALUE], predicate: Callable[[_KEY, _VALUE], bool]) -> Tuple[FSharpMap[_KEY, _VALUE], FSharpMap[_KEY, _VALUE]]:
    pattern_input: Tuple[Optional[MapTreeLeaf_2[_KEY, _VALUE]], Optional[MapTreeLeaf_2[_KEY, _VALUE]]] = MapTreeModule_partition(m.comparer, predicate, m.tree)
    return (FSharpMap__ctor(m.comparer, pattern_input[0]), FSharpMap__ctor(m.comparer, pattern_input[1]))


def FSharpMap__get_Count(m: FSharpMap[Any, Any]) -> int:
    return MapTreeModule_size(m.tree)


def FSharpMap__ContainsKey(m: FSharpMap[_KEY, Any], key: _KEY) -> bool:
    return MapTreeModule_mem(m.comparer, key, m.tree)


def FSharpMap__Remove(m: FSharpMap[_KEY, _VALUE], key: _KEY) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(m.comparer, MapTreeModule_remove(m.comparer, key, m.tree))


def FSharpMap__TryGetValue(_: FSharpMap[_KEY, _VALUE], key: _KEY, value: FSharpRef[_VALUE]) -> bool:
    match_value: Optional[_VALUE] = MapTreeModule_tryFind(_.comparer, key, _.tree)
    if match_value is None:
        return False

    else: 
        v: _VALUE = value_1(match_value)
        value.contents = v
        return True



def FSharpMap__get_Keys(__: FSharpMap[_KEY, Any]) -> ICollection[_KEY]:
    def mapping(kvp: Any, __: Any=__) -> _KEY:
        return kvp[0]

    return map_2(mapping, MapTreeModule_toArray(__.tree), None)


def FSharpMap__get_Values(__: FSharpMap[Any, _VALUE]) -> ICollection[_VALUE]:
    def mapping(kvp: Any, __: Any=__) -> _VALUE:
        return kvp[1]

    return map_2(mapping, MapTreeModule_toArray(__.tree), None)


def FSharpMap__TryFind(m: FSharpMap[_KEY, _VALUE], key: _KEY) -> Optional[_VALUE]:
    return MapTreeModule_tryFind(m.comparer, key, m.tree)


def FSharpMap__ToList(m: FSharpMap[_KEY, _VALUE]) -> FSharpList[Tuple[_KEY, _VALUE]]:
    return MapTreeModule_toList(m.tree)


def FSharpMap__ToArray(m: FSharpMap[_KEY, _VALUE]) -> Array[Any]:
    return MapTreeModule_toArray(m.tree)


def FSharpMap__ComputeHashCode(this: FSharpMap[Any, Any]) -> int:
    def combine_hash(x: int, y: int, this: Any=this) -> int:
        return ((x << 1) + y) + 631

    res: int = 0
    with get_enumerator(this) as enumerator:
        while enumerator.System_Collections_IEnumerator_MoveNext():
            active_pattern_result: Tuple[_KEY, _VALUE] = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
            res = combine_hash(res, structural_hash(active_pattern_result[0])) or 0
            res = combine_hash(res, structural_hash(active_pattern_result[1])) or 0
    return res


def is_empty(table: FSharpMap[Any, Any]) -> bool:
    return FSharpMap__get_IsEmpty(table)


def add(key: __A, value: __B, table: FSharpMap[__A, __B]) -> FSharpMap[__A, __B]:
    return FSharpMap__Add(table, key, value)


def change(key: __A, f: Callable[[Optional[__B]], Optional[__B]], table: FSharpMap[__A, __B]) -> FSharpMap[__A, __B]:
    return FSharpMap__Change(table, key, f)


def find(key: __A, table: FSharpMap[__A, __B]) -> __B:
    return FSharpMap__get_Item(table, key)


def try_find(key: __A, table: FSharpMap[__A, __B]) -> Optional[__B]:
    return FSharpMap__TryFind(table, key)


def remove(key: __A, table: FSharpMap[__A, __B]) -> FSharpMap[__A, __B]:
    return FSharpMap__Remove(table, key)


def contains_key(key: __A, table: FSharpMap[__A, Any]) -> bool:
    return FSharpMap__ContainsKey(table, key)


def iterate(action: Callable[[__A, __B], None], table: FSharpMap[__A, __B]) -> None:
    FSharpMap__Iterate(table, action)


def try_pick(chooser: Callable[[__A, __B], Optional[__C]], table: FSharpMap[__A, __B]) -> Optional[__C]:
    return FSharpMap__TryPick(table, chooser)


def pick(chooser: Callable[[__A, __B], Optional[__C]], table: FSharpMap[__A, __B]) -> __C:
    match_value: Optional[__C] = try_pick(chooser, table)
    if match_value is not None:
        return value_1(match_value)

    else: 
        raise Exception()



def exists(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> bool:
    return FSharpMap__Exists(table, predicate)


def filter(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> FSharpMap[__A, __B]:
    return FSharpMap__Filter(table, predicate)


def partition(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> Tuple[FSharpMap[__A, __B], FSharpMap[__A, __B]]:
    return FSharpMap__Partition(table, predicate)


def for_all(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> bool:
    return FSharpMap__ForAll(table, predicate)


def map(mapping: Callable[[__A, __B], __C], table: FSharpMap[__A, __B]) -> FSharpMap[__A, __C]:
    return FSharpMap__Map(table, mapping)


def fold(folder: Callable[[_STATE, _KEY, _T], _STATE], state: _STATE, table: FSharpMap[_KEY, _T]) -> _STATE:
    return MapTreeModule_fold(folder, state, FSharpMap__get_Tree(table))


def fold_back(folder: Callable[[_KEY, _T, _STATE], _STATE], table: FSharpMap[_KEY, _T], state: _STATE) -> _STATE:
    return MapTreeModule_foldBack(folder, FSharpMap__get_Tree(table), state)


def to_seq(table: FSharpMap[__A, __B]) -> IEnumerable_1[Tuple[__A, __B]]:
    def mapping(kvp: Any, table: Any=table) -> Tuple[__A, __B]:
        return (kvp[0], kvp[1])

    return map_1(mapping, table)


def find_key(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> __A:
    def chooser(kvp: Any, predicate: Any=predicate, table: Any=table) -> Optional[__A]:
        k: __A = kvp[0]
        if predicate(k, kvp[1]):
            return some(k)

        else: 
            return None


    return pick_1(chooser, table)


def try_find_key(predicate: Callable[[__A, __B], bool], table: FSharpMap[__A, __B]) -> Optional[__A]:
    def chooser(kvp: Any, predicate: Any=predicate, table: Any=table) -> Optional[__A]:
        k: __A = kvp[0]
        if predicate(k, kvp[1]):
            return some(k)

        else: 
            return None


    return try_pick_1(chooser, table)


def of_list(elements: FSharpList[Tuple[_KEY, _VALUE]], comparer: IComparer_1[_KEY]) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(comparer, MapTreeModule_ofSeq(comparer, elements))


def of_seq(elements: IEnumerable_1[Tuple[_T, __A]], comparer: IComparer_1[_T]) -> FSharpMap[_T, __A]:
    return FSharpMap__ctor(comparer, MapTreeModule_ofSeq(comparer, elements))


def of_array(elements: Array[Tuple[_KEY, _VALUE]], comparer: IComparer_1[_KEY]) -> FSharpMap[_KEY, _VALUE]:
    return FSharpMap__ctor(comparer, MapTreeModule_ofSeq(comparer, elements))


def to_list(table: FSharpMap[__A, __B]) -> FSharpList[Tuple[__A, __B]]:
    return FSharpMap__ToList(table)


def to_array(table: FSharpMap[__A, __B]) -> Array[Any]:
    return FSharpMap__ToArray(table)


def keys(table: FSharpMap[_K, Any]) -> ICollection[_K]:
    return FSharpMap__get_Keys(table)


def values(table: FSharpMap[Any, _V]) -> ICollection[_V]:
    return FSharpMap__get_Values(table)


def empty(comparer: IComparer_1[_KEY]) -> FSharpMap[_KEY, Any]:
    return FSharpMap_Empty(comparer)


def count(table: FSharpMap[Any, Any]) -> int:
    return FSharpMap__get_Count(table)


__all__ = ["MapTreeLeaf_2_reflection", "MapTreeLeaf_2__get_Key", "MapTreeLeaf_2__get_Value", "MapTreeNode_2_reflection", "MapTreeNode_2__get_Left", "MapTreeNode_2__get_Right", "MapTreeNode_2__get_Height", "MapTreeModule_empty", "MapTreeModule_sizeAux", "MapTreeModule_size", "MapTreeModule_mk", "MapTreeModule_rebalance", "MapTreeModule_add", "MapTreeModule_tryFind", "MapTreeModule_find", "MapTreeModule_partition1", "MapTreeModule_partitionAux", "MapTreeModule_partition", "MapTreeModule_filter1", "MapTreeModule_filterAux", "MapTreeModule_filter", "MapTreeModule_spliceOutSuccessor", "MapTreeModule_remove", "MapTreeModule_change", "MapTreeModule_mem", "MapTreeModule_iterOpt", "MapTreeModule_iter", "MapTreeModule_tryPickOpt", "MapTreeModule_tryPick", "MapTreeModule_existsOpt", "MapTreeModule_exists", "MapTreeModule_forallOpt", "MapTreeModule_forall", "MapTreeModule_map", "MapTreeModule_mapiOpt", "MapTreeModule_mapi", "MapTreeModule_foldBackOpt", "MapTreeModule_foldBack", "MapTreeModule_foldOpt", "MapTreeModule_fold", "MapTreeModule_foldSectionOpt", "MapTreeModule_foldSection", "MapTreeModule_toList", "MapTreeModule_copyToArray", "MapTreeModule_toArray", "MapTreeModule_ofList", "MapTreeModule_mkFromEnumerator", "MapTreeModule_ofArray", "MapTreeModule_ofSeq", "MapTreeModule_MapIterator_2_reflection", "MapTreeModule_collapseLHS", "MapTreeModule_mkIterator", "MapTreeModule_notStarted", "MapTreeModule_alreadyFinished", "MapTreeModule_current", "MapTreeModule_moveNext", "MapTreeModule_mkIEnumerator", "MapTreeModule_toSeq", "FSharpMap_reflection", "FSharpMap_Empty", "FSharpMap__get_Comparer", "FSharpMap__get_Tree", "FSharpMap__Add", "FSharpMap__Change", "FSharpMap__get_IsEmpty", "FSharpMap__get_Item", "FSharpMap__TryPick", "FSharpMap__Exists", "FSharpMap__Filter", "FSharpMap__ForAll", "FSharpMap__Fold", "FSharpMap__FoldSection", "FSharpMap__Iterate", "FSharpMap__MapRange", "FSharpMap__Map", "FSharpMap__Partition", "FSharpMap__get_Count", "FSharpMap__ContainsKey", "FSharpMap__Remove", "FSharpMap__TryGetValue", "FSharpMap__get_Keys", "FSharpMap__get_Values", "FSharpMap__TryFind", "FSharpMap__ToList", "FSharpMap__ToArray", "FSharpMap__ComputeHashCode", "is_empty", "add", "change", "find", "try_find", "remove", "contains_key", "iterate", "try_pick", "pick", "exists", "filter", "partition", "for_all", "map", "fold", "fold_back", "to_seq", "find_key", "try_find_key", "of_list", "of_seq", "of_array", "to_list", "to_array", "keys", "values", "empty", "count"]

