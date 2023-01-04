from __future__ import annotations
from dataclasses import dataclass
from typing import (TypeVar, Any, Optional, Generic, Tuple, Callable)
from .array import (fill, fold as fold_1)
from .list import (FSharpList, is_empty as is_empty_1, head, tail, of_array_with_tail, empty as empty_1, singleton as singleton_1, cons, fold as fold_2)
from .option import (value as value_1, some)
from .reflection import (TypeInfo, class_type, option_type, list_type, bool_type, record_type)
from .seq import (map as map_1, iterate as iterate_1, reduce, fold as fold_3, for_all as for_all_1, cache, exists as exists_1)
from .string import join
from .types import (Record, Array, to_string)
from .util import (IComparer_1, IEnumerator, IEnumerable_1, is_array_like, get_enumerator, ignore, to_iterator, structural_hash)
from .mutable_set import (HashSet, HashSet__ctor_Z6150332D, HashSet__get_Comparer)

_T = TypeVar("_T")

__A = TypeVar("__A")

__B = TypeVar("__B")

_A = TypeVar("_A")

_A_ = TypeVar("_A_")

_U = TypeVar("_U")

_STATE = TypeVar("_STATE")

def _expr231(gen0: TypeInfo) -> TypeInfo:
    return class_type("Set.SetTreeLeaf`1", [gen0], SetTreeLeaf_1)


class SetTreeLeaf_1(Generic[_T]):
    def __init__(self, k: Optional[Any]=None) -> None:
        self.k: _T = k


SetTreeLeaf_1_reflection = _expr231

def SetTreeLeaf_1__ctor_2B595(k: Optional[Any]=None) -> SetTreeLeaf_1[_T]:
    return SetTreeLeaf_1(k)


def SetTreeLeaf_1__get_Key(_: SetTreeLeaf_1[_T]) -> _T:
    return _.k


def _expr232(gen0: TypeInfo) -> TypeInfo:
    return class_type("Set.SetTreeNode`1", [gen0], SetTreeNode_1, SetTreeLeaf_1_reflection(gen0))


class SetTreeNode_1(SetTreeLeaf_1, Generic[_T]):
    def __init__(self, v: _T, left: Optional[SetTreeLeaf_1[_T]], right: Optional[SetTreeLeaf_1[_T]], h: int) -> None:
        super().__init__(v)
        self.left: Optional[SetTreeLeaf_1[_T]] = left
        self.right: Optional[SetTreeLeaf_1[_T]] = right
        self.h: int = h or 0


SetTreeNode_1_reflection = _expr232

def SetTreeNode_1__ctor_5F465FC9(v: _T, left: Optional[SetTreeLeaf_1[_T]], right: Optional[SetTreeLeaf_1[_T]], h: int) -> SetTreeNode_1[_T]:
    return SetTreeNode_1(v, left, right, h)


def SetTreeNode_1__get_Left(_: SetTreeNode_1[_T]) -> Optional[SetTreeLeaf_1[_T]]:
    return _.left


def SetTreeNode_1__get_Right(_: SetTreeNode_1[_T]) -> Optional[SetTreeLeaf_1[_T]]:
    return _.right


def SetTreeNode_1__get_Height(_: SetTreeNode_1[Any]) -> int:
    return _.h


def SetTreeModule_empty(__unit: None=None) -> Optional[SetTreeLeaf_1[_T]]:
    return None


def SetTreeModule_countAux(t_mut: Optional[SetTreeLeaf_1[_T]], acc_mut: int) -> int:
    while True:
        (t, acc) = (t_mut, acc_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                t_mut = SetTreeNode_1__get_Left(t2)
                acc_mut = SetTreeModule_countAux(SetTreeNode_1__get_Right(t2), acc + 1)
                continue

            else: 
                return acc + 1


        else: 
            return acc

        break


def SetTreeModule_count(s: Optional[SetTreeLeaf_1[__A]]=None) -> int:
    return SetTreeModule_countAux(s, 0)


def SetTreeModule_mk(l: Optional[SetTreeLeaf_1[_T]], k: _T, r: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    hl: int
    t: Optional[SetTreeLeaf_1[_T]] = l
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        hl = SetTreeNode_1__get_Height(t2) if isinstance(t2, SetTreeNode_1) else 1

    else: 
        hl = 0

    hr: int
    t_1: Optional[SetTreeLeaf_1[_T]] = r
    if t_1 is not None:
        t2_1: SetTreeLeaf_1[_T] = t_1
        hr = SetTreeNode_1__get_Height(t2_1) if isinstance(t2_1, SetTreeNode_1) else 1

    else: 
        hr = 0

    m: int = (hr if (hl < hr) else hl) or 0
    if m == 0:
        return SetTreeLeaf_1__ctor_2B595(k)

    else: 
        return SetTreeNode_1__ctor_5F465FC9(k, l, r, m + 1)



def SetTreeModule_rebalance(t1: Optional[SetTreeLeaf_1[_T]], v: _T, t2: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    t1h: int
    t: Optional[SetTreeLeaf_1[_T]] = t1
    if t is not None:
        t2_1: SetTreeLeaf_1[_T] = t
        t1h = SetTreeNode_1__get_Height(t2_1) if isinstance(t2_1, SetTreeNode_1) else 1

    else: 
        t1h = 0

    t2h: int
    t_1: Optional[SetTreeLeaf_1[_T]] = t2
    if t_1 is not None:
        t2_2: SetTreeLeaf_1[_T] = t_1
        t2h = SetTreeNode_1__get_Height(t2_2) if isinstance(t2_2, SetTreeNode_1) else 1

    else: 
        t2h = 0

    if t2h > (t1h + 2):
        match_value: SetTreeLeaf_1[_T] = value_1(t2)
        if isinstance(match_value, SetTreeNode_1):
            def _arrow233(__unit: None=None, t1: Any=t1, v: _T=v, t2: Any=t2) -> int:
                t_2: Optional[SetTreeLeaf_1[_T]] = SetTreeNode_1__get_Left(match_value)
                if t_2 is not None:
                    t2_3: SetTreeLeaf_1[_T] = t_2
                    return SetTreeNode_1__get_Height(t2_3) if isinstance(t2_3, SetTreeNode_1) else 1

                else: 
                    return 0


            if _arrow233() > (t1h + 1):
                match_value_1: SetTreeLeaf_1[_T] = value_1(SetTreeNode_1__get_Left(match_value))
                if isinstance(match_value_1, SetTreeNode_1):
                    return SetTreeModule_mk(SetTreeModule_mk(t1, v, SetTreeNode_1__get_Left(match_value_1)), SetTreeLeaf_1__get_Key(match_value_1), SetTreeModule_mk(SetTreeNode_1__get_Right(match_value_1), SetTreeLeaf_1__get_Key(match_value), SetTreeNode_1__get_Right(match_value)))

                else: 
                    raise Exception("internal error: Set.rebalance")


            else: 
                return SetTreeModule_mk(SetTreeModule_mk(t1, v, SetTreeNode_1__get_Left(match_value)), SetTreeLeaf_1__get_Key(match_value), SetTreeNode_1__get_Right(match_value))


        else: 
            raise Exception("internal error: Set.rebalance")


    elif t1h > (t2h + 2):
        match_value_2: SetTreeLeaf_1[_T] = value_1(t1)
        if isinstance(match_value_2, SetTreeNode_1):
            def _arrow234(__unit: None=None, t1: Any=t1, v: _T=v, t2: Any=t2) -> int:
                t_3: Optional[SetTreeLeaf_1[_T]] = SetTreeNode_1__get_Right(match_value_2)
                if t_3 is not None:
                    t2_4: SetTreeLeaf_1[_T] = t_3
                    return SetTreeNode_1__get_Height(t2_4) if isinstance(t2_4, SetTreeNode_1) else 1

                else: 
                    return 0


            if _arrow234() > (t2h + 1):
                match_value_3: SetTreeLeaf_1[_T] = value_1(SetTreeNode_1__get_Right(match_value_2))
                if isinstance(match_value_3, SetTreeNode_1):
                    return SetTreeModule_mk(SetTreeModule_mk(SetTreeNode_1__get_Left(match_value_2), SetTreeLeaf_1__get_Key(match_value_2), SetTreeNode_1__get_Left(match_value_3)), SetTreeLeaf_1__get_Key(match_value_3), SetTreeModule_mk(SetTreeNode_1__get_Right(match_value_3), v, t2))

                else: 
                    raise Exception("internal error: Set.rebalance")


            else: 
                return SetTreeModule_mk(SetTreeNode_1__get_Left(match_value_2), SetTreeLeaf_1__get_Key(match_value_2), SetTreeModule_mk(SetTreeNode_1__get_Right(match_value_2), v, t2))


        else: 
            raise Exception("internal error: Set.rebalance")


    else: 
        return SetTreeModule_mk(t1, v, t2)



def SetTreeModule_add(comparer: IComparer_1[_T], k: _T, t: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        c: int = comparer.Compare(k, SetTreeLeaf_1__get_Key(t2)) or 0
        if isinstance(t2, SetTreeNode_1):
            if c < 0:
                return SetTreeModule_rebalance(SetTreeModule_add(comparer, k, SetTreeNode_1__get_Left(t2)), SetTreeLeaf_1__get_Key(t2), SetTreeNode_1__get_Right(t2))

            elif c == 0:
                return t

            else: 
                return SetTreeModule_rebalance(SetTreeNode_1__get_Left(t2), SetTreeLeaf_1__get_Key(t2), SetTreeModule_add(comparer, k, SetTreeNode_1__get_Right(t2)))


        else: 
            c_1: int = comparer.Compare(k, SetTreeLeaf_1__get_Key(t2)) or 0
            if c_1 < 0:
                return SetTreeNode_1__ctor_5F465FC9(k, SetTreeModule_empty(), t, 2)

            elif c_1 == 0:
                return t

            else: 
                return SetTreeNode_1__ctor_5F465FC9(k, t, SetTreeModule_empty(), 2)



    else: 
        return SetTreeLeaf_1__ctor_2B595(k)



def SetTreeModule_balance(comparer: IComparer_1[_T], t1: Optional[SetTreeLeaf_1[_T]], k: _T, t2: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    if t1 is not None:
        t1_0027: SetTreeLeaf_1[_T] = t1
        if t2 is not None:
            t2_0027: SetTreeLeaf_1[_T] = t2
            if isinstance(t1_0027, SetTreeNode_1):
                if isinstance(t2_0027, SetTreeNode_1):
                    if (SetTreeNode_1__get_Height(t1_0027) + 2) < SetTreeNode_1__get_Height(t2_0027):
                        return SetTreeModule_rebalance(SetTreeModule_balance(comparer, t1, k, SetTreeNode_1__get_Left(t2_0027)), SetTreeLeaf_1__get_Key(t2_0027), SetTreeNode_1__get_Right(t2_0027))

                    elif (SetTreeNode_1__get_Height(t2_0027) + 2) < SetTreeNode_1__get_Height(t1_0027):
                        return SetTreeModule_rebalance(SetTreeNode_1__get_Left(t1_0027), SetTreeLeaf_1__get_Key(t1_0027), SetTreeModule_balance(comparer, SetTreeNode_1__get_Right(t1_0027), k, t2))

                    else: 
                        return SetTreeModule_mk(t1, k, t2)


                else: 
                    return SetTreeModule_add(comparer, k, SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2_0027), t1))


            else: 
                return SetTreeModule_add(comparer, k, SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t1_0027), t2))


        else: 
            return SetTreeModule_add(comparer, k, t1)


    else: 
        return SetTreeModule_add(comparer, k, t2)



def SetTreeModule_split(comparer: IComparer_1[_T], pivot: _T, t: Optional[SetTreeLeaf_1[_T]]=None) -> Tuple[Optional[SetTreeLeaf_1[_T]], bool, Optional[SetTreeLeaf_1[_T]]]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        if isinstance(t2, SetTreeNode_1):
            c: int = comparer.Compare(pivot, SetTreeLeaf_1__get_Key(t2)) or 0
            if c < 0:
                pattern_input: Tuple[Optional[SetTreeLeaf_1[_T]], bool, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_split(comparer, pivot, SetTreeNode_1__get_Left(t2))
                return (pattern_input[0], pattern_input[1], SetTreeModule_balance(comparer, pattern_input[2], SetTreeLeaf_1__get_Key(t2), SetTreeNode_1__get_Right(t2)))

            elif c == 0:
                return (SetTreeNode_1__get_Left(t2), True, SetTreeNode_1__get_Right(t2))

            else: 
                pattern_input_1: Tuple[Optional[SetTreeLeaf_1[_T]], bool, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_split(comparer, pivot, SetTreeNode_1__get_Right(t2))
                return (SetTreeModule_balance(comparer, SetTreeNode_1__get_Left(t2), SetTreeLeaf_1__get_Key(t2), pattern_input_1[0]), pattern_input_1[1], pattern_input_1[2])


        else: 
            c_1: int = comparer.Compare(SetTreeLeaf_1__get_Key(t2), pivot) or 0
            if c_1 < 0:
                return (t, False, SetTreeModule_empty())

            elif c_1 == 0:
                return (SetTreeModule_empty(), True, SetTreeModule_empty())

            else: 
                return (SetTreeModule_empty(), False, t)



    else: 
        return (SetTreeModule_empty(), False, SetTreeModule_empty())



def SetTreeModule_spliceOutSuccessor(t: Optional[SetTreeLeaf_1[_T]]=None) -> Tuple[_T, Optional[SetTreeLeaf_1[_T]]]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        if isinstance(t2, SetTreeNode_1):
            if SetTreeNode_1__get_Left(t2) is None:
                return (SetTreeLeaf_1__get_Key(t2), SetTreeNode_1__get_Right(t2))

            else: 
                pattern_input: Tuple[_T, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_spliceOutSuccessor(SetTreeNode_1__get_Left(t2))
                return (pattern_input[0], SetTreeModule_mk(pattern_input[1], SetTreeLeaf_1__get_Key(t2), SetTreeNode_1__get_Right(t2)))


        else: 
            return (SetTreeLeaf_1__get_Key(t2), SetTreeModule_empty())


    else: 
        raise Exception("internal error: Set.spliceOutSuccessor")



def SetTreeModule_remove(comparer: IComparer_1[_T], k: _T, t: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        c: int = comparer.Compare(k, SetTreeLeaf_1__get_Key(t2)) or 0
        if isinstance(t2, SetTreeNode_1):
            if c < 0:
                return SetTreeModule_rebalance(SetTreeModule_remove(comparer, k, SetTreeNode_1__get_Left(t2)), SetTreeLeaf_1__get_Key(t2), SetTreeNode_1__get_Right(t2))

            elif c == 0:
                if SetTreeNode_1__get_Left(t2) is None:
                    return SetTreeNode_1__get_Right(t2)

                elif SetTreeNode_1__get_Right(t2) is None:
                    return SetTreeNode_1__get_Left(t2)

                else: 
                    pattern_input: Tuple[_T, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_spliceOutSuccessor(SetTreeNode_1__get_Right(t2))
                    return SetTreeModule_mk(SetTreeNode_1__get_Left(t2), pattern_input[0], pattern_input[1])


            else: 
                return SetTreeModule_rebalance(SetTreeNode_1__get_Left(t2), SetTreeLeaf_1__get_Key(t2), SetTreeModule_remove(comparer, k, SetTreeNode_1__get_Right(t2)))


        elif c == 0:
            return SetTreeModule_empty()

        else: 
            return t


    else: 
        return t



def SetTreeModule_mem(comparer_mut: IComparer_1[_T], k_mut: _T, t_mut: Optional[SetTreeLeaf_1[_T]]) -> bool:
    while True:
        (comparer, k, t) = (comparer_mut, k_mut, t_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            c: int = comparer.Compare(k, SetTreeLeaf_1__get_Key(t2)) or 0
            if isinstance(t2, SetTreeNode_1):
                if c < 0:
                    comparer_mut = comparer
                    k_mut = k
                    t_mut = SetTreeNode_1__get_Left(t2)
                    continue

                elif c == 0:
                    return True

                else: 
                    comparer_mut = comparer
                    k_mut = k
                    t_mut = SetTreeNode_1__get_Right(t2)
                    continue


            else: 
                return c == 0


        else: 
            return False

        break


def SetTreeModule_iter(f_mut: Callable[[_T], None], t_mut: Optional[SetTreeLeaf_1[_T]]) -> None:
    while True:
        (f, t) = (f_mut, t_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                SetTreeModule_iter(f, SetTreeNode_1__get_Left(t2))
                f(SetTreeLeaf_1__get_Key(t2))
                f_mut = f
                t_mut = SetTreeNode_1__get_Right(t2)
                continue

            else: 
                f(SetTreeLeaf_1__get_Key(t2))


        break


def SetTreeModule_foldBackOpt(f_mut: Any, t_mut: Optional[SetTreeLeaf_1[_T]], x_mut: __A) -> __A:
    while True:
        (f, t, x) = (f_mut, t_mut, x_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                f_mut = f
                t_mut = SetTreeNode_1__get_Left(t2)
                x_mut = f(SetTreeLeaf_1__get_Key(t2), SetTreeModule_foldBackOpt(f, SetTreeNode_1__get_Right(t2), x))
                continue

            else: 
                return f(SetTreeLeaf_1__get_Key(t2), x)


        else: 
            return x

        break


def SetTreeModule_foldBack(f: Callable[[__A, __B], __B], m: Optional[SetTreeLeaf_1[__A]], x: __B) -> __B:
    return SetTreeModule_foldBackOpt(f, m, x)


def SetTreeModule_foldOpt(f_mut: Any, x_mut: __A, t_mut: Optional[SetTreeLeaf_1[_T]]) -> __A:
    while True:
        (f, x, t) = (f_mut, x_mut, t_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                f_mut = f
                x_mut = f(SetTreeModule_foldOpt(f, x, SetTreeNode_1__get_Left(t2)), SetTreeLeaf_1__get_Key(t2))
                t_mut = SetTreeNode_1__get_Right(t2)
                continue

            else: 
                return f(x, SetTreeLeaf_1__get_Key(t2))


        else: 
            return x

        break


def SetTreeModule_fold(f: Callable[[__A, __B], __A], x: __A, m: Optional[SetTreeLeaf_1[__B]]=None) -> __A:
    return SetTreeModule_foldOpt(f, x, m)


def SetTreeModule_forall(f_mut: Callable[[_T], bool], t_mut: Optional[SetTreeLeaf_1[_T]]) -> bool:
    while True:
        (f, t) = (f_mut, t_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                if SetTreeModule_forall(f, SetTreeNode_1__get_Left(t2)) if f(SetTreeLeaf_1__get_Key(t2)) else False:
                    f_mut = f
                    t_mut = SetTreeNode_1__get_Right(t2)
                    continue

                else: 
                    return False


            else: 
                return f(SetTreeLeaf_1__get_Key(t2))


        else: 
            return True

        break


def SetTreeModule_exists(f_mut: Callable[[_T], bool], t_mut: Optional[SetTreeLeaf_1[_T]]) -> bool:
    while True:
        (f, t) = (f_mut, t_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                if True if f(SetTreeLeaf_1__get_Key(t2)) else SetTreeModule_exists(f, SetTreeNode_1__get_Left(t2)):
                    return True

                else: 
                    f_mut = f
                    t_mut = SetTreeNode_1__get_Right(t2)
                    continue


            else: 
                return f(SetTreeLeaf_1__get_Key(t2))


        else: 
            return False

        break


def SetTreeModule_subset(comparer: IComparer_1[__A], a: Optional[SetTreeLeaf_1[__A]]=None, b: Optional[SetTreeLeaf_1[__A]]=None) -> bool:
    def _arrow235(x: Optional[__A]=None, comparer: IComparer_1[Any]=comparer, a: Any=a, b: Any=b) -> bool:
        return SetTreeModule_mem(comparer, x, b)

    return SetTreeModule_forall(_arrow235, a)


def SetTreeModule_properSubset(comparer: IComparer_1[__A], a: Optional[SetTreeLeaf_1[__A]]=None, b: Optional[SetTreeLeaf_1[__A]]=None) -> bool:
    def _arrow236(x: Optional[__A]=None, comparer: IComparer_1[Any]=comparer, a: Any=a, b: Any=b) -> bool:
        return SetTreeModule_mem(comparer, x, b)

    if SetTreeModule_forall(_arrow236, a):
        def _arrow237(x_1: Optional[__A]=None, comparer: IComparer_1[Any]=comparer, a: Any=a, b: Any=b) -> bool:
            return not SetTreeModule_mem(comparer, x_1, a)

        return SetTreeModule_exists(_arrow237, b)

    else: 
        return False



def SetTreeModule_filterAux(comparer_mut: IComparer_1[_T], f_mut: Callable[[_T], bool], t_mut: Optional[SetTreeLeaf_1[_T]], acc_mut: Optional[SetTreeLeaf_1[_T]]) -> Optional[SetTreeLeaf_1[_T]]:
    while True:
        (comparer, f, t, acc) = (comparer_mut, f_mut, t_mut, acc_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                acc_1: Optional[SetTreeLeaf_1[_T]] = SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2), acc) if f(SetTreeLeaf_1__get_Key(t2)) else acc
                comparer_mut = comparer
                f_mut = f
                t_mut = SetTreeNode_1__get_Left(t2)
                acc_mut = SetTreeModule_filterAux(comparer, f, SetTreeNode_1__get_Right(t2), acc_1)
                continue

            elif f(SetTreeLeaf_1__get_Key(t2)):
                return SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2), acc)

            else: 
                return acc


        else: 
            return acc

        break


def SetTreeModule_filter(comparer: IComparer_1[__A], f: Callable[[__A], bool], s: Optional[SetTreeLeaf_1[__A]]=None) -> Optional[SetTreeLeaf_1[__A]]:
    return SetTreeModule_filterAux(comparer, f, s, SetTreeModule_empty())


def SetTreeModule_diffAux(comparer_mut: IComparer_1[_T], t_mut: Optional[SetTreeLeaf_1[_T]], acc_mut: Optional[SetTreeLeaf_1[_T]]) -> Optional[SetTreeLeaf_1[_T]]:
    while True:
        (comparer, t, acc) = (comparer_mut, t_mut, acc_mut)
        if acc is None:
            return acc

        elif t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                comparer_mut = comparer
                t_mut = SetTreeNode_1__get_Left(t2)
                acc_mut = SetTreeModule_diffAux(comparer, SetTreeNode_1__get_Right(t2), SetTreeModule_remove(comparer, SetTreeLeaf_1__get_Key(t2), acc))
                continue

            else: 
                return SetTreeModule_remove(comparer, SetTreeLeaf_1__get_Key(t2), acc)


        else: 
            return acc

        break


def SetTreeModule_diff(comparer: IComparer_1[__A], a: Optional[SetTreeLeaf_1[__A]]=None, b: Optional[SetTreeLeaf_1[__A]]=None) -> Optional[SetTreeLeaf_1[__A]]:
    return SetTreeModule_diffAux(comparer, b, a)


def SetTreeModule_union(comparer: IComparer_1[_T], t1: Optional[SetTreeLeaf_1[_T]]=None, t2: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[SetTreeLeaf_1[_T]]:
    if t1 is not None:
        t1_0027: SetTreeLeaf_1[_T] = t1
        if t2 is not None:
            t2_0027: SetTreeLeaf_1[_T] = t2
            if isinstance(t1_0027, SetTreeNode_1):
                if isinstance(t2_0027, SetTreeNode_1):
                    if SetTreeNode_1__get_Height(t1_0027) > SetTreeNode_1__get_Height(t2_0027):
                        pattern_input: Tuple[Optional[SetTreeLeaf_1[_T]], bool, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_split(comparer, SetTreeLeaf_1__get_Key(t1_0027), t2)
                        return SetTreeModule_balance(comparer, SetTreeModule_union(comparer, SetTreeNode_1__get_Left(t1_0027), pattern_input[0]), SetTreeLeaf_1__get_Key(t1_0027), SetTreeModule_union(comparer, SetTreeNode_1__get_Right(t1_0027), pattern_input[2]))

                    else: 
                        pattern_input_1: Tuple[Optional[SetTreeLeaf_1[_T]], bool, Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_split(comparer, SetTreeLeaf_1__get_Key(t2_0027), t1)
                        return SetTreeModule_balance(comparer, SetTreeModule_union(comparer, SetTreeNode_1__get_Left(t2_0027), pattern_input_1[0]), SetTreeLeaf_1__get_Key(t2_0027), SetTreeModule_union(comparer, SetTreeNode_1__get_Right(t2_0027), pattern_input_1[2]))


                else: 
                    return SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2_0027), t1)


            else: 
                return SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t1_0027), t2)


        else: 
            return t1


    else: 
        return t2



def SetTreeModule_intersectionAux(comparer_mut: IComparer_1[_T], b_mut: Optional[SetTreeLeaf_1[_T]], t_mut: Optional[SetTreeLeaf_1[_T]], acc_mut: Optional[SetTreeLeaf_1[_T]]) -> Optional[SetTreeLeaf_1[_T]]:
    while True:
        (comparer, b, t, acc) = (comparer_mut, b_mut, t_mut, acc_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                acc_1: Optional[SetTreeLeaf_1[_T]] = SetTreeModule_intersectionAux(comparer, b, SetTreeNode_1__get_Right(t2), acc)
                acc_2: Optional[SetTreeLeaf_1[_T]] = SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2), acc_1) if SetTreeModule_mem(comparer, SetTreeLeaf_1__get_Key(t2), b) else acc_1
                comparer_mut = comparer
                b_mut = b
                t_mut = SetTreeNode_1__get_Left(t2)
                acc_mut = acc_2
                continue

            elif SetTreeModule_mem(comparer, SetTreeLeaf_1__get_Key(t2), b):
                return SetTreeModule_add(comparer, SetTreeLeaf_1__get_Key(t2), acc)

            else: 
                return acc


        else: 
            return acc

        break


def SetTreeModule_intersection(comparer: IComparer_1[__A], a: Optional[SetTreeLeaf_1[__A]]=None, b: Optional[SetTreeLeaf_1[__A]]=None) -> Optional[SetTreeLeaf_1[__A]]:
    return SetTreeModule_intersectionAux(comparer, b, a, SetTreeModule_empty())


def SetTreeModule_partition1(comparer: IComparer_1[__A], f: Callable[[__A], bool], k: __A, acc1: Optional[SetTreeLeaf_1[__A]]=None, acc2: Optional[SetTreeLeaf_1[__A]]=None) -> Tuple[Optional[SetTreeLeaf_1[__A]], Optional[SetTreeLeaf_1[__A]]]:
    if f(k):
        return (SetTreeModule_add(comparer, k, acc1), acc2)

    else: 
        return (acc1, SetTreeModule_add(comparer, k, acc2))



def SetTreeModule_partitionAux(comparer_mut: IComparer_1[_T], f_mut: Callable[[_T], bool], t_mut: Optional[SetTreeLeaf_1[_T]], acc__mut: Optional[SetTreeLeaf_1[_T]], acc__1_mut: Optional[SetTreeLeaf_1[_T]]) -> Tuple[Optional[SetTreeLeaf_1[_T]], Optional[SetTreeLeaf_1[_T]]]:
    while True:
        (comparer, f, t, acc_, acc__1) = (comparer_mut, f_mut, t_mut, acc__mut, acc__1_mut)
        acc: Tuple[Optional[SetTreeLeaf_1[_T]], Optional[SetTreeLeaf_1[_T]]] = (acc_, acc__1)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                acc_1: Tuple[Optional[SetTreeLeaf_1[_T]], Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_partitionAux(comparer, f, SetTreeNode_1__get_Right(t2), acc[0], acc[1])
                acc_4: Tuple[Optional[SetTreeLeaf_1[_T]], Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_partition1(comparer, f, SetTreeLeaf_1__get_Key(t2), acc_1[0], acc_1[1])
                comparer_mut = comparer
                f_mut = f
                t_mut = SetTreeNode_1__get_Left(t2)
                acc__mut = acc_4[0]
                acc__1_mut = acc_4[1]
                continue

            else: 
                return SetTreeModule_partition1(comparer, f, SetTreeLeaf_1__get_Key(t2), acc[0], acc[1])


        else: 
            return acc

        break


def SetTreeModule_partition(comparer: IComparer_1[__A], f: Callable[[__A], bool], s: Optional[SetTreeLeaf_1[__A]]=None) -> Tuple[Optional[SetTreeLeaf_1[__A]], Optional[SetTreeLeaf_1[__A]]]:
    return SetTreeModule_partitionAux(comparer, f, s, SetTreeModule_empty(), SetTreeModule_empty())


def SetTreeModule_minimumElementAux(t_mut: Optional[SetTreeLeaf_1[_T]], n_mut: _T) -> _T:
    while True:
        (t, n) = (t_mut, n_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                t_mut = SetTreeNode_1__get_Left(t2)
                n_mut = SetTreeLeaf_1__get_Key(t2)
                continue

            else: 
                return SetTreeLeaf_1__get_Key(t2)


        else: 
            return n

        break


def SetTreeModule_minimumElementOpt(t: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[_T]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        if isinstance(t2, SetTreeNode_1):
            return some(SetTreeModule_minimumElementAux(SetTreeNode_1__get_Left(t2), SetTreeLeaf_1__get_Key(t2)))

        else: 
            return some(SetTreeLeaf_1__get_Key(t2))


    else: 
        return None



def SetTreeModule_maximumElementAux(t_mut: Optional[SetTreeLeaf_1[_T]], n_mut: _T) -> _T:
    while True:
        (t, n) = (t_mut, n_mut)
        if t is not None:
            t2: SetTreeLeaf_1[_T] = t
            if isinstance(t2, SetTreeNode_1):
                t_mut = SetTreeNode_1__get_Right(t2)
                n_mut = SetTreeLeaf_1__get_Key(t2)
                continue

            else: 
                return SetTreeLeaf_1__get_Key(t2)


        else: 
            return n

        break


def SetTreeModule_maximumElementOpt(t: Optional[SetTreeLeaf_1[_T]]=None) -> Optional[_T]:
    if t is not None:
        t2: SetTreeLeaf_1[_T] = t
        if isinstance(t2, SetTreeNode_1):
            return some(SetTreeModule_maximumElementAux(SetTreeNode_1__get_Right(t2), SetTreeLeaf_1__get_Key(t2)))

        else: 
            return some(SetTreeLeaf_1__get_Key(t2))


    else: 
        return None



def SetTreeModule_minimumElement(s: Optional[SetTreeLeaf_1[__A]]=None) -> __A:
    match_value: Optional[__A] = SetTreeModule_minimumElementOpt(s)
    if match_value is None:
        raise Exception("Set contains no elements")

    else: 
        return value_1(match_value)



def SetTreeModule_maximumElement(s: Optional[SetTreeLeaf_1[__A]]=None) -> __A:
    match_value: Optional[__A] = SetTreeModule_maximumElementOpt(s)
    if match_value is None:
        raise Exception("Set contains no elements")

    else: 
        return value_1(match_value)



def _expr238(gen0: TypeInfo) -> TypeInfo:
    return record_type("Set.SetTreeModule.SetIterator`1", [gen0], SetTreeModule_SetIterator_1, lambda: [("stack", list_type(option_type(SetTreeLeaf_1_reflection(gen0)))), ("started", bool_type)])


@dataclass(eq = False, repr = False)
class SetTreeModule_SetIterator_1(Record, Generic[_T]):
    stack: FSharpList[Optional[SetTreeLeaf_1[_T]]]
    started: bool

SetTreeModule_SetIterator_1_reflection = _expr238

def SetTreeModule_collapseLHS(stack_mut: FSharpList[Optional[SetTreeLeaf_1[_T]]]) -> FSharpList[Optional[SetTreeLeaf_1[_T]]]:
    while True:
        (stack,) = (stack_mut,)
        if not is_empty_1(stack):
            x: Optional[SetTreeLeaf_1[_T]] = head(stack)
            rest: FSharpList[Optional[SetTreeLeaf_1[_T]]] = tail(stack)
            if x is not None:
                x2: SetTreeLeaf_1[_T] = x
                if isinstance(x2, SetTreeNode_1):
                    stack_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2)), SetTreeNode_1__get_Right(x2)], rest)
                    continue

                else: 
                    return stack


            else: 
                stack_mut = rest
                continue


        else: 
            return empty_1()

        break


def SetTreeModule_mkIterator(s: Optional[SetTreeLeaf_1[__A]]=None) -> SetTreeModule_SetIterator_1[__A]:
    return SetTreeModule_SetIterator_1(SetTreeModule_collapseLHS(singleton_1(s)), False)


def SetTreeModule_notStarted(__unit: None=None) -> Any:
    raise Exception("Enumeration not started")


def SetTreeModule_alreadyFinished(__unit: None=None) -> Any:
    raise Exception("Enumeration already started")


def SetTreeModule_current(i: SetTreeModule_SetIterator_1[__A]) -> __A:
    if i.started:
        match_value: FSharpList[Optional[SetTreeLeaf_1[__A]]] = i.stack
        if is_empty_1(match_value):
            return SetTreeModule_alreadyFinished()

        elif head(match_value) is not None:
            t: SetTreeLeaf_1[__A] = head(match_value)
            return SetTreeLeaf_1__get_Key(t)

        else: 
            raise Exception("Please report error: Set iterator, unexpected stack for current")


    else: 
        return SetTreeModule_notStarted()



def SetTreeModule_moveNext(i: SetTreeModule_SetIterator_1[Any]) -> bool:
    if i.started:
        match_value: FSharpList[Optional[SetTreeLeaf_1[_T]]] = i.stack
        if not is_empty_1(match_value):
            if head(match_value) is not None:
                t: SetTreeLeaf_1[_T] = head(match_value)
                if isinstance(t, SetTreeNode_1):
                    raise Exception("Please report error: Set iterator, unexpected stack for moveNext")

                else: 
                    i.stack = SetTreeModule_collapseLHS(tail(match_value))
                    return not is_empty_1(i.stack)


            else: 
                raise Exception("Please report error: Set iterator, unexpected stack for moveNext")


        else: 
            return False


    else: 
        i.started = True
        return not is_empty_1(i.stack)



def SetTreeModule_mkIEnumerator(s: Optional[SetTreeLeaf_1[_A]]=None) -> IEnumerator[_A]:
    i: SetTreeModule_SetIterator_1[_A] = SetTreeModule_mkIterator(s)
    class ObjectExpr239(IEnumerator[_A_]):
        def System_Collections_Generic_IEnumerator_1_get_Current(self, __unit: None=None, s: Any=s) -> _A_:
            return SetTreeModule_current(i)

        def System_Collections_IEnumerator_get_Current(self, __unit: None=None, s: Any=s) -> Any:
            return SetTreeModule_current(i)

        def System_Collections_IEnumerator_MoveNext(self, __unit: None=None, s: Any=s) -> bool:
            return SetTreeModule_moveNext(i)

        def System_Collections_IEnumerator_Reset(self, __unit: None=None, s: Any=s) -> None:
            nonlocal i
            i = SetTreeModule_mkIterator(s)

        def Dispose(self, __unit: None=None, s: Any=s) -> None:
            pass

    return ObjectExpr239()


def SetTreeModule_compareStacks(comparer_mut: IComparer_1[_T], l1_mut: FSharpList[Optional[SetTreeLeaf_1[_T]]], l2_mut: FSharpList[Optional[SetTreeLeaf_1[_T]]]) -> int:
    while True:
        (comparer, l1, l2) = (comparer_mut, l1_mut, l2_mut)
        if not is_empty_1(l1):
            if not is_empty_1(l2):
                if head(l2) is not None:
                    if head(l1) is not None:
                        x1_3: SetTreeLeaf_1[_T] = head(l1)
                        x2_3: SetTreeLeaf_1[_T] = head(l2)
                        if isinstance(x1_3, SetTreeNode_1):
                            if SetTreeNode_1__get_Left(x1_3) is None:
                                if isinstance(x2_3, SetTreeNode_1):
                                    if SetTreeNode_1__get_Left(x2_3) is None:
                                        c: int = comparer.Compare(SetTreeLeaf_1__get_Key(x1_3), SetTreeLeaf_1__get_Key(x2_3)) or 0
                                        if c != 0:
                                            return c

                                        else: 
                                            comparer_mut = comparer
                                            l1_mut = cons(SetTreeNode_1__get_Right(x1_3), tail(l1))
                                            l2_mut = cons(SetTreeNode_1__get_Right(x2_3), tail(l2))
                                            continue


                                    else: 
                                        (pattern_matching_result, t1_6, x1_4, t2_6, x2_4) = (None, None, None, None, None)
                                        if not is_empty_1(l1):
                                            if head(l1) is not None:
                                                pattern_matching_result = 0
                                                t1_6 = tail(l1)
                                                x1_4 = head(l1)

                                            elif not is_empty_1(l2):
                                                if head(l2) is not None:
                                                    pattern_matching_result = 1
                                                    t2_6 = tail(l2)
                                                    x2_4 = head(l2)

                                                else: 
                                                    pattern_matching_result = 2


                                            else: 
                                                pattern_matching_result = 2


                                        elif not is_empty_1(l2):
                                            if head(l2) is not None:
                                                pattern_matching_result = 1
                                                t2_6 = tail(l2)
                                                x2_4 = head(l2)

                                            else: 
                                                pattern_matching_result = 2


                                        else: 
                                            pattern_matching_result = 2

                                        if pattern_matching_result == 0:
                                            if isinstance(x1_4, SetTreeNode_1):
                                                comparer_mut = comparer
                                                l1_mut = of_array_with_tail([SetTreeNode_1__get_Left(x1_4), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x1_4), SetTreeModule_empty(), SetTreeNode_1__get_Right(x1_4), 0)], t1_6)
                                                l2_mut = l2
                                                continue

                                            else: 
                                                comparer_mut = comparer
                                                l1_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x1_4))], t1_6)
                                                l2_mut = l2
                                                continue


                                        elif pattern_matching_result == 1:
                                            if isinstance(x2_4, SetTreeNode_1):
                                                comparer_mut = comparer
                                                l1_mut = l1
                                                l2_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2_4), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x2_4), SetTreeModule_empty(), SetTreeNode_1__get_Right(x2_4), 0)], t2_6)
                                                continue

                                            else: 
                                                comparer_mut = comparer
                                                l1_mut = l1
                                                l2_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2_4))], t2_6)
                                                continue


                                        elif pattern_matching_result == 2:
                                            raise Exception("unexpected state in SetTree.compareStacks")



                                else: 
                                    c_1: int = comparer.Compare(SetTreeLeaf_1__get_Key(x1_3), SetTreeLeaf_1__get_Key(x2_3)) or 0
                                    if c_1 != 0:
                                        return c_1

                                    else: 
                                        comparer_mut = comparer
                                        l1_mut = cons(SetTreeNode_1__get_Right(x1_3), tail(l1))
                                        l2_mut = cons(SetTreeModule_empty(), tail(l2))
                                        continue



                            else: 
                                (pattern_matching_result_1, t1_7, x1_5, t2_7, x2_5) = (None, None, None, None, None)
                                if not is_empty_1(l1):
                                    if head(l1) is not None:
                                        pattern_matching_result_1 = 0
                                        t1_7 = tail(l1)
                                        x1_5 = head(l1)

                                    elif not is_empty_1(l2):
                                        if head(l2) is not None:
                                            pattern_matching_result_1 = 1
                                            t2_7 = tail(l2)
                                            x2_5 = head(l2)

                                        else: 
                                            pattern_matching_result_1 = 2


                                    else: 
                                        pattern_matching_result_1 = 2


                                elif not is_empty_1(l2):
                                    if head(l2) is not None:
                                        pattern_matching_result_1 = 1
                                        t2_7 = tail(l2)
                                        x2_5 = head(l2)

                                    else: 
                                        pattern_matching_result_1 = 2


                                else: 
                                    pattern_matching_result_1 = 2

                                if pattern_matching_result_1 == 0:
                                    if isinstance(x1_5, SetTreeNode_1):
                                        comparer_mut = comparer
                                        l1_mut = of_array_with_tail([SetTreeNode_1__get_Left(x1_5), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x1_5), SetTreeModule_empty(), SetTreeNode_1__get_Right(x1_5), 0)], t1_7)
                                        l2_mut = l2
                                        continue

                                    else: 
                                        comparer_mut = comparer
                                        l1_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x1_5))], t1_7)
                                        l2_mut = l2
                                        continue


                                elif pattern_matching_result_1 == 1:
                                    if isinstance(x2_5, SetTreeNode_1):
                                        comparer_mut = comparer
                                        l1_mut = l1
                                        l2_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2_5), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x2_5), SetTreeModule_empty(), SetTreeNode_1__get_Right(x2_5), 0)], t2_7)
                                        continue

                                    else: 
                                        comparer_mut = comparer
                                        l1_mut = l1
                                        l2_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2_5))], t2_7)
                                        continue


                                elif pattern_matching_result_1 == 2:
                                    raise Exception("unexpected state in SetTree.compareStacks")



                        elif isinstance(x2_3, SetTreeNode_1):
                            if SetTreeNode_1__get_Left(x2_3) is None:
                                c_2: int = comparer.Compare(SetTreeLeaf_1__get_Key(x1_3), SetTreeLeaf_1__get_Key(x2_3)) or 0
                                if c_2 != 0:
                                    return c_2

                                else: 
                                    comparer_mut = comparer
                                    l1_mut = cons(SetTreeModule_empty(), tail(l1))
                                    l2_mut = cons(SetTreeNode_1__get_Right(x2_3), tail(l2))
                                    continue


                            else: 
                                (pattern_matching_result_2, t1_8, x1_6, t2_8, x2_6) = (None, None, None, None, None)
                                if not is_empty_1(l1):
                                    if head(l1) is not None:
                                        pattern_matching_result_2 = 0
                                        t1_8 = tail(l1)
                                        x1_6 = head(l1)

                                    elif not is_empty_1(l2):
                                        if head(l2) is not None:
                                            pattern_matching_result_2 = 1
                                            t2_8 = tail(l2)
                                            x2_6 = head(l2)

                                        else: 
                                            pattern_matching_result_2 = 2


                                    else: 
                                        pattern_matching_result_2 = 2


                                elif not is_empty_1(l2):
                                    if head(l2) is not None:
                                        pattern_matching_result_2 = 1
                                        t2_8 = tail(l2)
                                        x2_6 = head(l2)

                                    else: 
                                        pattern_matching_result_2 = 2


                                else: 
                                    pattern_matching_result_2 = 2

                                if pattern_matching_result_2 == 0:
                                    if isinstance(x1_6, SetTreeNode_1):
                                        comparer_mut = comparer
                                        l1_mut = of_array_with_tail([SetTreeNode_1__get_Left(x1_6), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x1_6), SetTreeModule_empty(), SetTreeNode_1__get_Right(x1_6), 0)], t1_8)
                                        l2_mut = l2
                                        continue

                                    else: 
                                        comparer_mut = comparer
                                        l1_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x1_6))], t1_8)
                                        l2_mut = l2
                                        continue


                                elif pattern_matching_result_2 == 1:
                                    if isinstance(x2_6, SetTreeNode_1):
                                        comparer_mut = comparer
                                        l1_mut = l1
                                        l2_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2_6), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x2_6), SetTreeModule_empty(), SetTreeNode_1__get_Right(x2_6), 0)], t2_8)
                                        continue

                                    else: 
                                        comparer_mut = comparer
                                        l1_mut = l1
                                        l2_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2_6))], t2_8)
                                        continue


                                elif pattern_matching_result_2 == 2:
                                    raise Exception("unexpected state in SetTree.compareStacks")



                        else: 
                            c_3: int = comparer.Compare(SetTreeLeaf_1__get_Key(x1_3), SetTreeLeaf_1__get_Key(x2_3)) or 0
                            if c_3 != 0:
                                return c_3

                            else: 
                                comparer_mut = comparer
                                l1_mut = tail(l1)
                                l2_mut = tail(l2)
                                continue



                    else: 
                        x2: SetTreeLeaf_1[_T] = head(l2)
                        (pattern_matching_result_3, t1_2, x1, t2_2, x2_1) = (None, None, None, None, None)
                        if not is_empty_1(l1):
                            if head(l1) is not None:
                                pattern_matching_result_3 = 0
                                t1_2 = tail(l1)
                                x1 = head(l1)

                            elif not is_empty_1(l2):
                                if head(l2) is not None:
                                    pattern_matching_result_3 = 1
                                    t2_2 = tail(l2)
                                    x2_1 = head(l2)

                                else: 
                                    pattern_matching_result_3 = 2


                            else: 
                                pattern_matching_result_3 = 2


                        elif not is_empty_1(l2):
                            if head(l2) is not None:
                                pattern_matching_result_3 = 1
                                t2_2 = tail(l2)
                                x2_1 = head(l2)

                            else: 
                                pattern_matching_result_3 = 2


                        else: 
                            pattern_matching_result_3 = 2

                        if pattern_matching_result_3 == 0:
                            if isinstance(x1, SetTreeNode_1):
                                comparer_mut = comparer
                                l1_mut = of_array_with_tail([SetTreeNode_1__get_Left(x1), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x1), SetTreeModule_empty(), SetTreeNode_1__get_Right(x1), 0)], t1_2)
                                l2_mut = l2
                                continue

                            else: 
                                comparer_mut = comparer
                                l1_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x1))], t1_2)
                                l2_mut = l2
                                continue


                        elif pattern_matching_result_3 == 1:
                            if isinstance(x2_1, SetTreeNode_1):
                                comparer_mut = comparer
                                l1_mut = l1
                                l2_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2_1), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x2_1), SetTreeModule_empty(), SetTreeNode_1__get_Right(x2_1), 0)], t2_2)
                                continue

                            else: 
                                comparer_mut = comparer
                                l1_mut = l1
                                l2_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2_1))], t2_2)
                                continue


                        elif pattern_matching_result_3 == 2:
                            raise Exception("unexpected state in SetTree.compareStacks")



                elif head(l1) is not None:
                    x1_1: SetTreeLeaf_1[_T] = head(l1)
                    (pattern_matching_result_4, t1_4, x1_2, t2_4, x2_2) = (None, None, None, None, None)
                    if not is_empty_1(l1):
                        if head(l1) is not None:
                            pattern_matching_result_4 = 0
                            t1_4 = tail(l1)
                            x1_2 = head(l1)

                        elif not is_empty_1(l2):
                            if head(l2) is not None:
                                pattern_matching_result_4 = 1
                                t2_4 = tail(l2)
                                x2_2 = head(l2)

                            else: 
                                pattern_matching_result_4 = 2


                        else: 
                            pattern_matching_result_4 = 2


                    elif not is_empty_1(l2):
                        if head(l2) is not None:
                            pattern_matching_result_4 = 1
                            t2_4 = tail(l2)
                            x2_2 = head(l2)

                        else: 
                            pattern_matching_result_4 = 2


                    else: 
                        pattern_matching_result_4 = 2

                    if pattern_matching_result_4 == 0:
                        if isinstance(x1_2, SetTreeNode_1):
                            comparer_mut = comparer
                            l1_mut = of_array_with_tail([SetTreeNode_1__get_Left(x1_2), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x1_2), SetTreeModule_empty(), SetTreeNode_1__get_Right(x1_2), 0)], t1_4)
                            l2_mut = l2
                            continue

                        else: 
                            comparer_mut = comparer
                            l1_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x1_2))], t1_4)
                            l2_mut = l2
                            continue


                    elif pattern_matching_result_4 == 1:
                        if isinstance(x2_2, SetTreeNode_1):
                            comparer_mut = comparer
                            l1_mut = l1
                            l2_mut = of_array_with_tail([SetTreeNode_1__get_Left(x2_2), SetTreeNode_1__ctor_5F465FC9(SetTreeLeaf_1__get_Key(x2_2), SetTreeModule_empty(), SetTreeNode_1__get_Right(x2_2), 0)], t2_4)
                            continue

                        else: 
                            comparer_mut = comparer
                            l1_mut = l1
                            l2_mut = of_array_with_tail([SetTreeModule_empty(), SetTreeLeaf_1__ctor_2B595(SetTreeLeaf_1__get_Key(x2_2))], t2_4)
                            continue


                    elif pattern_matching_result_4 == 2:
                        raise Exception("unexpected state in SetTree.compareStacks")


                else: 
                    comparer_mut = comparer
                    l1_mut = tail(l1)
                    l2_mut = tail(l2)
                    continue


            else: 
                return 1


        elif is_empty_1(l2):
            return 0

        else: 
            return -1

        break


def SetTreeModule_compare(comparer: IComparer_1[_T], t1: Optional[SetTreeLeaf_1[_T]]=None, t2: Optional[SetTreeLeaf_1[_T]]=None) -> int:
    if t1 is None:
        if t2 is None:
            return 0

        else: 
            return -1


    elif t2 is None:
        return 1

    else: 
        return SetTreeModule_compareStacks(comparer, singleton_1(t1), singleton_1(t2))



def SetTreeModule_choose(s: Optional[SetTreeLeaf_1[__A]]=None) -> __A:
    return SetTreeModule_minimumElement(s)


def SetTreeModule_toList(t: Optional[SetTreeLeaf_1[_T]]=None) -> FSharpList[_T]:
    def loop(t_0027_mut: Optional[SetTreeLeaf_1[_T]], acc_mut: FSharpList[_T], t: Any=t) -> FSharpList[_T]:
        while True:
            (t_0027, acc) = (t_0027_mut, acc_mut)
            if t_0027 is not None:
                t2: SetTreeLeaf_1[_T] = t_0027
                if isinstance(t2, SetTreeNode_1):
                    t_0027_mut = SetTreeNode_1__get_Left(t2)
                    acc_mut = cons(SetTreeLeaf_1__get_Key(t2), loop(SetTreeNode_1__get_Right(t2), acc))
                    continue

                else: 
                    return cons(SetTreeLeaf_1__get_Key(t2), acc)


            else: 
                return acc

            break

    return loop(t, empty_1())


def SetTreeModule_copyToArray(s: Optional[SetTreeLeaf_1[__A]], arr: Array[__A], i: int) -> None:
    j: int = i or 0
    def _arrow240(x: Optional[__A]=None, s: Any=s, arr: Array[Any]=arr, i: int=i) -> None:
        nonlocal j
        arr[j] = x
        j = (j + 1) or 0

    SetTreeModule_iter(_arrow240, s)


def SetTreeModule_toArray(s: Optional[SetTreeLeaf_1[__A]]=None) -> Array[__A]:
    n: int = SetTreeModule_count(s) or 0
    res: Array[__A] = fill([0] * n, 0, n, None)
    SetTreeModule_copyToArray(s, res, 0)
    return res


def SetTreeModule_mkFromEnumerator(comparer_mut: IComparer_1[__A], acc_mut: Optional[SetTreeLeaf_1[__A]], e_mut: IEnumerator[__A]) -> Optional[SetTreeLeaf_1[__A]]:
    while True:
        (comparer, acc, e) = (comparer_mut, acc_mut, e_mut)
        if e.System_Collections_IEnumerator_MoveNext():
            comparer_mut = comparer
            acc_mut = SetTreeModule_add(comparer, e.System_Collections_Generic_IEnumerator_1_get_Current(), acc)
            e_mut = e
            continue

        else: 
            return acc

        break


def SetTreeModule_ofArray(comparer: IComparer_1[__A], l: Array[__A]) -> Optional[SetTreeLeaf_1[__A]]:
    def _arrow241(acc: Optional[SetTreeLeaf_1[__A]], k: __A, comparer: IComparer_1[Any]=comparer, l: Array[Any]=l) -> Optional[SetTreeLeaf_1[__A]]:
        return SetTreeModule_add(comparer, k, acc)

    return fold_1(_arrow241, SetTreeModule_empty(), l)


def SetTreeModule_ofList(comparer: IComparer_1[__A], l: FSharpList[__A]) -> Optional[SetTreeLeaf_1[__A]]:
    def _arrow242(acc: Optional[SetTreeLeaf_1[__A]], k: __A, comparer: IComparer_1[Any]=comparer, l: FSharpList[Any]=l) -> Optional[SetTreeLeaf_1[__A]]:
        return SetTreeModule_add(comparer, k, acc)

    return fold_2(_arrow242, SetTreeModule_empty(), l)


def SetTreeModule_ofSeq(comparer: IComparer_1[_T], c: IEnumerable_1[_T]) -> Optional[SetTreeLeaf_1[_T]]:
    if is_array_like(c):
        return SetTreeModule_ofArray(comparer, c)

    elif isinstance(c, FSharpList):
        return SetTreeModule_ofList(comparer, c)

    else: 
        with get_enumerator(c) as ie:
            return SetTreeModule_mkFromEnumerator(comparer, SetTreeModule_empty(), ie)



def _expr244(gen0: TypeInfo) -> TypeInfo:
    return class_type("Set.FSharpSet", [gen0], FSharpSet)


class FSharpSet(Generic[_T]):
    def __init__(self, comparer: IComparer_1[_T], tree: Optional[SetTreeLeaf_1[_T]]=None) -> None:
        self.comparer: IComparer_1[_T] = comparer
        self.tree: Optional[SetTreeLeaf_1[_T]] = tree

    def GetHashCode(self, __unit: None=None) -> int:
        this: FSharpSet[_T] = self
        return FSharpSet__ComputeHashCode(this)

    def __eq__(self, that: Any=None) -> bool:
        this: FSharpSet[_T] = self
        return (SetTreeModule_compare(FSharpSet__get_Comparer(this), FSharpSet__get_Tree(this), FSharpSet__get_Tree(that)) == 0) if isinstance(that, FSharpSet) else False

    def __str__(self, __unit: None=None) -> str:
        this: FSharpSet[_T] = self
        def _arrow243(x: Optional[_T]=None) -> str:
            copy_of_struct: _T = x
            return to_string(copy_of_struct)

        return ("set [" + join("; ", map_1(_arrow243, this))) + "]"

    @property
    def Symbol_toStringTag(self, __unit: None=None) -> str:
        return "FSharpSet"

    def to_json(self, __unit: None=None) -> Any:
        this: FSharpSet[_T] = self
        return list(this)

    def __cmp__(self, that: Any=None) -> int:
        s: FSharpSet[_T] = self
        return SetTreeModule_compare(FSharpSet__get_Comparer(s), FSharpSet__get_Tree(s), FSharpSet__get_Tree(that))

    def System_Collections_Generic_ICollection_1_Add2B595(self, x: Optional[_T]=None) -> None:
        ignore(x)
        raise Exception("ReadOnlyCollection")

    def System_Collections_Generic_ICollection_1_Clear(self, __unit: None=None) -> None:
        raise Exception("ReadOnlyCollection")

    def System_Collections_Generic_ICollection_1_Remove2B595(self, x: Optional[_T]=None) -> bool:
        ignore(x)
        raise Exception("ReadOnlyCollection")

    def System_Collections_Generic_ICollection_1_Contains2B595(self, x: Optional[_T]=None) -> bool:
        s: FSharpSet[_T] = self
        return SetTreeModule_mem(FSharpSet__get_Comparer(s), x, FSharpSet__get_Tree(s))

    def System_Collections_Generic_ICollection_1_CopyToZ3B4C077E(self, arr: Array[_T], i: int) -> None:
        s: FSharpSet[_T] = self
        SetTreeModule_copyToArray(FSharpSet__get_Tree(s), arr, i)

    def System_Collections_Generic_ICollection_1_get_IsReadOnly(self, __unit: None=None) -> bool:
        return True

    def __len__(self, __unit: None=None) -> int:
        s: FSharpSet[_T] = self
        return FSharpSet__get_Count(s)

    def __len__(self, __unit: None=None) -> int:
        s: FSharpSet[_T] = self
        return FSharpSet__get_Count(s)

    def GetEnumerator(self, __unit: None=None) -> IEnumerator[_T]:
        s: FSharpSet[_T] = self
        return SetTreeModule_mkIEnumerator(FSharpSet__get_Tree(s))

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        s: FSharpSet[_T] = self
        return SetTreeModule_mkIEnumerator(FSharpSet__get_Tree(s))

    @property
    def size(self, __unit: None=None) -> int:
        s: FSharpSet[_T] = self
        return FSharpSet__get_Count(s)

    def add(self, k: Optional[_T]=None) -> Set_1[_T]:
        s: FSharpSet[_T] = self
        raise Exception("Set cannot be mutated")
        return s

    def clear(self, __unit: None=None) -> None:
        raise Exception("Set cannot be mutated")

    def delete(self, k: Optional[_T]=None) -> bool:
        raise Exception("Set cannot be mutated")
        return False

    def __contains__(self, k: Optional[_T]=None) -> bool:
        s: FSharpSet[_T] = self
        return FSharpSet__Contains(s, k)

    def keys(self, __unit: None=None) -> IEnumerable_1[_T]:
        s: FSharpSet[_T] = self
        def mapping(x: Optional[_T]=None) -> Optional[_T]:
            return x

        return map_1(mapping, s)

    def values(self, __unit: None=None) -> IEnumerable_1[_T]:
        s: FSharpSet[_T] = self
        def mapping(x: Optional[_T]=None) -> Optional[_T]:
            return x

        return map_1(mapping, s)

    def entries(self, __unit: None=None) -> IEnumerable_1[Tuple[_T, _T]]:
        s: FSharpSet[_T] = self
        def mapping(v: Optional[_T]=None) -> Tuple[_T, _T]:
            return (v, v)

        return map_1(mapping, s)

    def for_each(self, f: Callable[[_T, _T, Set_1[_T]], None], this_arg: Optional[Any]=None) -> None:
        s: FSharpSet[_T] = self
        def action(x: Optional[_T]=None) -> None:
            f(x, x, s)

        iterate_1(action, s)


FSharpSet_reflection = _expr244

def FSharpSet__ctor(comparer: IComparer_1[_T], tree: Optional[SetTreeLeaf_1[_T]]=None) -> FSharpSet[_T]:
    return FSharpSet(comparer, tree)


def FSharpSet__get_Comparer(set_1: FSharpSet[_T]) -> IComparer_1[_T]:
    return set_1.comparer


def FSharpSet__get_Tree(set_1: FSharpSet[_T]) -> Optional[SetTreeLeaf_1[_T]]:
    return set_1.tree


def FSharpSet_Empty(comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet__ctor(comparer, SetTreeModule_empty())


def FSharpSet__Add(s: FSharpSet[_T], value: _T) -> FSharpSet[_T]:
    return FSharpSet__ctor(FSharpSet__get_Comparer(s), SetTreeModule_add(FSharpSet__get_Comparer(s), value, FSharpSet__get_Tree(s)))


def FSharpSet__Remove(s: FSharpSet[_T], value: _T) -> FSharpSet[_T]:
    return FSharpSet__ctor(FSharpSet__get_Comparer(s), SetTreeModule_remove(FSharpSet__get_Comparer(s), value, FSharpSet__get_Tree(s)))


def FSharpSet__get_Count(s: FSharpSet[Any]) -> int:
    return SetTreeModule_count(FSharpSet__get_Tree(s))


def FSharpSet__Contains(s: FSharpSet[_T], value: _T) -> bool:
    return SetTreeModule_mem(FSharpSet__get_Comparer(s), value, FSharpSet__get_Tree(s))


def FSharpSet__Iterate(s: FSharpSet[_T], x: Callable[[_T], None]) -> None:
    SetTreeModule_iter(x, FSharpSet__get_Tree(s))


def FSharpSet__Fold(s: FSharpSet[_T], f: Callable[[_T, __A], __A], z: __A) -> __A:
    f_1: Any = f
    def _arrow245(x: __A, z_1: _T, s: Any=s, f: Any=f, z: __A=z) -> __A:
        return f_1(z_1, x)

    return SetTreeModule_fold(_arrow245, z, FSharpSet__get_Tree(s))


def FSharpSet__get_IsEmpty(s: FSharpSet[Any]) -> bool:
    return FSharpSet__get_Tree(s) is None


def FSharpSet__Partition(s: FSharpSet[_T], f: Callable[[_T], bool]) -> Tuple[FSharpSet[_T], FSharpSet[_T]]:
    if FSharpSet__get_Tree(s) is None:
        return (s, s)

    else: 
        pattern_input: Tuple[Optional[SetTreeLeaf_1[_T]], Optional[SetTreeLeaf_1[_T]]] = SetTreeModule_partition(FSharpSet__get_Comparer(s), f, FSharpSet__get_Tree(s))
        return (FSharpSet__ctor(FSharpSet__get_Comparer(s), pattern_input[0]), FSharpSet__ctor(FSharpSet__get_Comparer(s), pattern_input[1]))



def FSharpSet__Filter(s: FSharpSet[_T], f: Callable[[_T], bool]) -> FSharpSet[_T]:
    if FSharpSet__get_Tree(s) is None:
        return s

    else: 
        return FSharpSet__ctor(FSharpSet__get_Comparer(s), SetTreeModule_filter(FSharpSet__get_Comparer(s), f, FSharpSet__get_Tree(s)))



def FSharpSet__Map(s: FSharpSet[_T], f: Callable[[_T], _U], comparer: IComparer_1[_U]) -> FSharpSet[_U]:
    def _arrow246(acc: Optional[SetTreeLeaf_1[_U]], k: _T, s: Any=s, f: Any=f, comparer: IComparer_1[Any]=comparer) -> Optional[SetTreeLeaf_1[_U]]:
        return SetTreeModule_add(comparer, f(k), acc)

    return FSharpSet__ctor(comparer, SetTreeModule_fold(_arrow246, SetTreeModule_empty(), FSharpSet__get_Tree(s)))


def FSharpSet__Exists(s: FSharpSet[_T], f: Callable[[_T], bool]) -> bool:
    return SetTreeModule_exists(f, FSharpSet__get_Tree(s))


def FSharpSet__ForAll(s: FSharpSet[_T], f: Callable[[_T], bool]) -> bool:
    return SetTreeModule_forall(f, FSharpSet__get_Tree(s))


def FSharpSet_op_Subtraction(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> FSharpSet[_T]:
    if FSharpSet__get_Tree(set1) is None:
        return set1

    elif FSharpSet__get_Tree(set2) is None:
        return set1

    else: 
        return FSharpSet__ctor(FSharpSet__get_Comparer(set1), SetTreeModule_diff(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set1), FSharpSet__get_Tree(set2)))



def FSharpSet_op_Addition(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> FSharpSet[_T]:
    if FSharpSet__get_Tree(set2) is None:
        return set1

    elif FSharpSet__get_Tree(set1) is None:
        return set2

    else: 
        return FSharpSet__ctor(FSharpSet__get_Comparer(set1), SetTreeModule_union(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set1), FSharpSet__get_Tree(set2)))



def FSharpSet_Intersection(a: FSharpSet[_T], b: FSharpSet[_T]) -> FSharpSet[_T]:
    if FSharpSet__get_Tree(b) is None:
        return b

    elif FSharpSet__get_Tree(a) is None:
        return a

    else: 
        return FSharpSet__ctor(FSharpSet__get_Comparer(a), SetTreeModule_intersection(FSharpSet__get_Comparer(a), FSharpSet__get_Tree(a), FSharpSet__get_Tree(b)))



def FSharpSet_IntersectionMany(sets: IEnumerable_1[FSharpSet[_T]]) -> FSharpSet[_T]:
    def _arrow247(s1: FSharpSet[_T], s2: FSharpSet[_T], sets: Any=sets) -> FSharpSet[_T]:
        return FSharpSet_Intersection(s1, s2)

    return reduce(_arrow247, sets)


def FSharpSet_Equality(a: FSharpSet[_T], b: FSharpSet[_T]) -> bool:
    return SetTreeModule_compare(FSharpSet__get_Comparer(a), FSharpSet__get_Tree(a), FSharpSet__get_Tree(b)) == 0


def FSharpSet_Compare(a: FSharpSet[_T], b: FSharpSet[_T]) -> int:
    return SetTreeModule_compare(FSharpSet__get_Comparer(a), FSharpSet__get_Tree(a), FSharpSet__get_Tree(b))


def FSharpSet__get_Choose(x: FSharpSet[_T]) -> _T:
    return SetTreeModule_choose(FSharpSet__get_Tree(x))


def FSharpSet__get_MinimumElement(x: FSharpSet[_T]) -> _T:
    return SetTreeModule_minimumElement(FSharpSet__get_Tree(x))


def FSharpSet__get_MaximumElement(x: FSharpSet[_T]) -> _T:
    return SetTreeModule_maximumElement(FSharpSet__get_Tree(x))


def FSharpSet__IsSubsetOf(x: FSharpSet[_T], other_set: FSharpSet[_T]) -> bool:
    return SetTreeModule_subset(FSharpSet__get_Comparer(x), FSharpSet__get_Tree(x), FSharpSet__get_Tree(other_set))


def FSharpSet__IsSupersetOf(x: FSharpSet[_T], other_set: FSharpSet[_T]) -> bool:
    return SetTreeModule_subset(FSharpSet__get_Comparer(x), FSharpSet__get_Tree(other_set), FSharpSet__get_Tree(x))


def FSharpSet__IsProperSubsetOf(x: FSharpSet[_T], other_set: FSharpSet[_T]) -> bool:
    return SetTreeModule_properSubset(FSharpSet__get_Comparer(x), FSharpSet__get_Tree(x), FSharpSet__get_Tree(other_set))


def FSharpSet__IsProperSupersetOf(x: FSharpSet[_T], other_set: FSharpSet[_T]) -> bool:
    return SetTreeModule_properSubset(FSharpSet__get_Comparer(x), FSharpSet__get_Tree(other_set), FSharpSet__get_Tree(x))


def FSharpSet__ToList(x: FSharpSet[_T]) -> FSharpList[_T]:
    return SetTreeModule_toList(FSharpSet__get_Tree(x))


def FSharpSet__ToArray(x: FSharpSet[_T]) -> Array[_T]:
    return SetTreeModule_toArray(FSharpSet__get_Tree(x))


def FSharpSet__ComputeHashCode(this: FSharpSet[Any]) -> int:
    def combine_hash(x: int, y: int, this: Any=this) -> int:
        return ((x << 1) + y) + 631

    res: int = 0
    with get_enumerator(this) as enumerator:
        while enumerator.System_Collections_IEnumerator_MoveNext():
            x_1: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
            res = combine_hash(res, structural_hash(x_1)) or 0
    return abs(res)


def is_empty(set_1: FSharpSet[Any]) -> bool:
    return FSharpSet__get_IsEmpty(set_1)


def contains(element: _T, set_1: FSharpSet[_T]) -> bool:
    return FSharpSet__Contains(set_1, element)


def add(value: _T, set_1: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet__Add(set_1, value)


def singleton(value: _T, comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet__Add(FSharpSet_Empty(comparer), value)


def remove(value: _T, set_1: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet__Remove(set_1, value)


def union(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet_op_Addition(set1, set2)


def union_many(sets: IEnumerable_1[FSharpSet[_T]], comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    def _arrow248(s1: FSharpSet[_T], s2: FSharpSet[_T], sets: Any=sets, comparer: IComparer_1[Any]=comparer) -> FSharpSet[_T]:
        return FSharpSet_op_Addition(s1, s2)

    return fold_3(_arrow248, FSharpSet_Empty(comparer), sets)


def intersect(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet_Intersection(set1, set2)


def intersect_many(sets: IEnumerable_1[FSharpSet[_T]]) -> FSharpSet[_T]:
    return FSharpSet_IntersectionMany(sets)


def iterate(action: Callable[[_T], None], set_1: FSharpSet[_T]) -> None:
    FSharpSet__Iterate(set_1, action)


def empty(comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet_Empty(comparer)


def for_all(predicate: Callable[[_T], bool], set_1: FSharpSet[_T]) -> bool:
    return FSharpSet__ForAll(set_1, predicate)


def exists(predicate: Callable[[_T], bool], set_1: FSharpSet[_T]) -> bool:
    return FSharpSet__Exists(set_1, predicate)


def filter(predicate: Callable[[_T], bool], set_1: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet__Filter(set_1, predicate)


def partition(predicate: Callable[[_T], bool], set_1: FSharpSet[_T]) -> Tuple[FSharpSet[_T], FSharpSet[_T]]:
    return FSharpSet__Partition(set_1, predicate)


def fold(folder: Callable[[_STATE, _T], _STATE], state: _STATE, set_1: FSharpSet[_T]) -> _STATE:
    return SetTreeModule_fold(folder, state, FSharpSet__get_Tree(set_1))


def fold_back(folder: Callable[[_T, _STATE], _STATE], set_1: FSharpSet[_T], state: _STATE) -> _STATE:
    return SetTreeModule_foldBack(folder, FSharpSet__get_Tree(set_1), state)


def map(mapping: Callable[[_T], _U], set_1: FSharpSet[_T], comparer: IComparer_1[_U]) -> FSharpSet[_U]:
    return FSharpSet__Map(set_1, mapping, comparer)


def count(set_1: FSharpSet[Any]) -> int:
    return FSharpSet__get_Count(set_1)


def of_list(elements: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet__ctor(comparer, SetTreeModule_ofSeq(comparer, elements))


def of_array(array: Array[_T], comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet__ctor(comparer, SetTreeModule_ofArray(comparer, array))


def to_list(set_1: FSharpSet[_T]) -> FSharpList[_T]:
    return FSharpSet__ToList(set_1)


def to_array(set_1: FSharpSet[_T]) -> Array[_T]:
    return FSharpSet__ToArray(set_1)


def to_seq(set_1: FSharpSet[_T]) -> IEnumerable_1[_T]:
    def mapping(x: Optional[_T]=None, set_1: Any=set_1) -> Optional[_T]:
        return x

    return map_1(mapping, set_1)


def of_seq(elements: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> FSharpSet[_T]:
    return FSharpSet__ctor(comparer, SetTreeModule_ofSeq(comparer, elements))


def difference(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> FSharpSet[_T]:
    return FSharpSet_op_Subtraction(set1, set2)


def is_subset(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> bool:
    return SetTreeModule_subset(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set1), FSharpSet__get_Tree(set2))


def is_superset(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> bool:
    return SetTreeModule_subset(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set2), FSharpSet__get_Tree(set1))


def is_proper_subset(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> bool:
    return SetTreeModule_properSubset(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set1), FSharpSet__get_Tree(set2))


def is_proper_superset(set1: FSharpSet[_T], set2: FSharpSet[_T]) -> bool:
    return SetTreeModule_properSubset(FSharpSet__get_Comparer(set1), FSharpSet__get_Tree(set2), FSharpSet__get_Tree(set1))


def min_element(set_1: FSharpSet[_T]) -> _T:
    return FSharpSet__get_MinimumElement(set_1)


def max_element(set_1: FSharpSet[_T]) -> _T:
    return FSharpSet__get_MaximumElement(set_1)


def union_with(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> Set_1[_T]:
    def folder(acc: Set_1[_T], x: _T, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> Set_1[_T]:
        return acc.add(x)

    return fold_3(folder, s1, s2)


def new_mutable_set_with(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> Set_1[_T]:
    if isinstance(s1, HashSet):
        return HashSet__ctor_Z6150332D(s2, HashSet__get_Comparer(s1))

    else: 
        return Set(s2)



def intersect_with(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> None:
    s2_1: Set_1[_T] = new_mutable_set_with(s1, s2)
    def action(x: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> None:
        if not (x in s2_1):
            ignore(s1.delete(x))


    iterate_1(action, s1.values())


def except_with(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> None:
    def action(x: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> None:
        ignore(s1.delete(x))

    iterate_1(action, s2)


def is_subset_of(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> bool:
    s2_1: Set_1[_T] = new_mutable_set_with(s1, s2)
    def predicate(arg: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> bool:
        return arg in s2_1

    return for_all_1(predicate, s1.values())


def is_superset_of(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> bool:
    def predicate(arg: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> bool:
        return arg in s1

    return for_all_1(predicate, s2)


def is_proper_subset_of(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> bool:
    s2_1: Set_1[_T] = new_mutable_set_with(s1, s2)
    if s2_1.size > s1.size:
        def predicate(arg: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> bool:
            return arg in s2_1

        return for_all_1(predicate, s1.values())

    else: 
        return False



def is_proper_superset_of(s1: Set_1[_T], s2: IEnumerable_1[_T]) -> bool:
    s2_1: IEnumerable_1[_T] = cache(s2)
    def predicate(arg_1: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> bool:
        return not (arg_1 in s1)

    if exists_1(predicate, s2_1):
        def predicate_1(arg_2: Optional[_T]=None, s1: Any=s1, s2: IEnumerable_1[Any]=s2) -> bool:
            return arg_2 in s1

        return for_all_1(predicate_1, s2_1)

    else: 
        return False



__all__ = ["SetTreeLeaf_1_reflection", "SetTreeLeaf_1__get_Key", "SetTreeNode_1_reflection", "SetTreeNode_1__get_Left", "SetTreeNode_1__get_Right", "SetTreeNode_1__get_Height", "SetTreeModule_empty", "SetTreeModule_countAux", "SetTreeModule_count", "SetTreeModule_mk", "SetTreeModule_rebalance", "SetTreeModule_add", "SetTreeModule_balance", "SetTreeModule_split", "SetTreeModule_spliceOutSuccessor", "SetTreeModule_remove", "SetTreeModule_mem", "SetTreeModule_iter", "SetTreeModule_foldBackOpt", "SetTreeModule_foldBack", "SetTreeModule_foldOpt", "SetTreeModule_fold", "SetTreeModule_forall", "SetTreeModule_exists", "SetTreeModule_subset", "SetTreeModule_properSubset", "SetTreeModule_filterAux", "SetTreeModule_filter", "SetTreeModule_diffAux", "SetTreeModule_diff", "SetTreeModule_union", "SetTreeModule_intersectionAux", "SetTreeModule_intersection", "SetTreeModule_partition1", "SetTreeModule_partitionAux", "SetTreeModule_partition", "SetTreeModule_minimumElementAux", "SetTreeModule_minimumElementOpt", "SetTreeModule_maximumElementAux", "SetTreeModule_maximumElementOpt", "SetTreeModule_minimumElement", "SetTreeModule_maximumElement", "SetTreeModule_SetIterator_1_reflection", "SetTreeModule_collapseLHS", "SetTreeModule_mkIterator", "SetTreeModule_notStarted", "SetTreeModule_alreadyFinished", "SetTreeModule_current", "SetTreeModule_moveNext", "SetTreeModule_mkIEnumerator", "SetTreeModule_compareStacks", "SetTreeModule_compare", "SetTreeModule_choose", "SetTreeModule_toList", "SetTreeModule_copyToArray", "SetTreeModule_toArray", "SetTreeModule_mkFromEnumerator", "SetTreeModule_ofArray", "SetTreeModule_ofList", "SetTreeModule_ofSeq", "FSharpSet_reflection", "FSharpSet__get_Comparer", "FSharpSet__get_Tree", "FSharpSet_Empty", "FSharpSet__Add", "FSharpSet__Remove", "FSharpSet__get_Count", "FSharpSet__Contains", "FSharpSet__Iterate", "FSharpSet__Fold", "FSharpSet__get_IsEmpty", "FSharpSet__Partition", "FSharpSet__Filter", "FSharpSet__Map", "FSharpSet__Exists", "FSharpSet__ForAll", "FSharpSet_op_Subtraction", "FSharpSet_op_Addition", "FSharpSet_Intersection", "FSharpSet_IntersectionMany", "FSharpSet_Equality", "FSharpSet_Compare", "FSharpSet__get_Choose", "FSharpSet__get_MinimumElement", "FSharpSet__get_MaximumElement", "FSharpSet__IsSubsetOf", "FSharpSet__IsSupersetOf", "FSharpSet__IsProperSubsetOf", "FSharpSet__IsProperSupersetOf", "FSharpSet__ToList", "FSharpSet__ToArray", "FSharpSet__ComputeHashCode", "is_empty", "contains", "add", "singleton", "remove", "union", "union_many", "intersect", "intersect_many", "iterate", "empty", "for_all", "exists", "filter", "partition", "fold", "fold_back", "map", "count", "of_list", "of_array", "to_list", "to_array", "to_seq", "of_seq", "difference", "is_subset", "is_superset", "is_proper_subset", "is_proper_superset", "min_element", "max_element", "union_with", "new_mutable_set_with", "intersect_with", "except_with", "is_subset_of", "is_superset_of", "is_proper_subset_of", "is_proper_superset_of"]

