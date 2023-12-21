from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Generic, TypeVar

from .array_ import chunk_by_size as chunk_by_size_1
from .array_ import fill, sort_in_place_with
from .array_ import fold_back as fold_back_1
from .array_ import fold_back2 as fold_back2_1
from .array_ import iterate as iterate_1
from .array_ import map as map_1
from .array_ import pairwise as pairwise_1
from .array_ import permute as permute_1
from .array_ import scan_back as scan_back_1
from .array_ import split_into as split_into_1
from .array_ import transpose as transpose_1
from .array_ import try_find_back as try_find_back_1
from .array_ import try_find_index_back as try_find_index_back_1
from .array_ import windowed as windowed_1
from .global_ import (
    IGenericAdder_1,
    IGenericAverager_1,
    SR_differentLengths,
    SR_indexOutOfBounds,
    SR_inputMustBeNonNegative,
    SR_inputSequenceEmpty,
    SR_inputSequenceTooLong,
    SR_inputWasEmpty,
    SR_keyNotFoundAlt,
    SR_notEnoughElements,
)
from .option import default_arg, some
from .option import value as value_1
from .reflection import TypeInfo, class_type, option_type, record_type
from .string_ import join
from .types import Array, Record
from .util import (
    IComparer_1,
    IDisposable,
    IEnumerable_1,
    IEnumerator,
    IEqualityComparer_1,
    compare,
    equals,
    get_enumerator,
    ignore,
    is_array_like,
    structural_hash,
    to_iterator,
)


_T = TypeVar("_T")

__A = TypeVar("__A")

_STATE = TypeVar("_STATE")

_T1 = TypeVar("_T1")

_T2 = TypeVar("_T2")

__B = TypeVar("__B")

_U = TypeVar("_U")

_T3 = TypeVar("_T3")

_RESULT = TypeVar("_RESULT")

__C = TypeVar("__C")


def _expr68(gen0: TypeInfo) -> TypeInfo:
    return record_type(
        "ListModule.FSharpList",
        [gen0],
        FSharpList,
        lambda: [("head", gen0), ("tail", option_type(FSharpList_reflection(gen0)))],
    )


@dataclass(eq=False, repr=False, slots=True)
class FSharpList(Record, Generic[_T]):
    head: _T
    tail: FSharpList[_T] | None

    def __str__(self, __unit: None = None) -> str:
        xs: FSharpList[_T] = self
        return ("[" + join("; ", xs)) + "]"

    def __eq__(self, other: Any = None) -> bool:
        xs: FSharpList[_T] = self
        if xs is other:
            return True

        else:

            def loop(xs_1_mut: FSharpList[_T], ys_1_mut: FSharpList[_T]) -> bool:
                while True:
                    (xs_1, ys_1) = (xs_1_mut, ys_1_mut)
                    matchValue: FSharpList[_T] | None = xs_1.tail
                    matchValue_1: FSharpList[_T] | None = ys_1.tail
                    if matchValue is not None:
                        if matchValue_1 is not None:
                            xt: FSharpList[_T] = matchValue
                            yt: FSharpList[_T] = matchValue_1
                            if equals(xs_1.head, ys_1.head):
                                xs_1_mut = xt
                                ys_1_mut = yt
                                continue

                            else:
                                return False

                        else:
                            return False

                    elif matchValue_1 is not None:
                        return False

                    else:
                        return True

                    break

            return loop(xs, other)

    def GetHashCode(self, __unit: None = None) -> int:
        xs: FSharpList[_T] = self

        def loop(i_mut: int, h_mut: int, xs_1_mut: FSharpList[_T]) -> int:
            while True:
                (i, h, xs_1) = (i_mut, h_mut, xs_1_mut)
                match_value: FSharpList[_T] | None = xs_1.tail
                if match_value is not None:
                    t: FSharpList[_T] = match_value
                    if i > 18:
                        return h

                    else:
                        i_mut = i + 1
                        h_mut = ((h << 1) + structural_hash(xs_1.head)) + (631 * i)
                        xs_1_mut = t
                        continue

                else:
                    return h

                break

        return loop(0, 0, xs)

    def to_json(self, __unit: None = None) -> Any:
        this: FSharpList[_T] = self
        return list(this)

    def __cmp__(self, other: Any = None) -> int:
        xs: FSharpList[_T] = self

        def loop(xs_1_mut: FSharpList[_T], ys_1_mut: FSharpList[_T]) -> int:
            while True:
                (xs_1, ys_1) = (xs_1_mut, ys_1_mut)
                matchValue: FSharpList[_T] | None = xs_1.tail
                matchValue_1: FSharpList[_T] | None = ys_1.tail
                if matchValue is not None:
                    if matchValue_1 is not None:
                        xt: FSharpList[_T] = matchValue
                        yt: FSharpList[_T] = matchValue_1
                        c: int = compare(xs_1.head, ys_1.head) or 0
                        if c == 0:
                            xs_1_mut = xt
                            ys_1_mut = yt
                            continue

                        else:
                            return c

                    else:
                        return 1

                elif matchValue_1 is not None:
                    return -1

                else:
                    return 0

                break

        return loop(xs, other)

    def GetEnumerator(self, __unit: None = None) -> IEnumerator[_T]:
        xs: FSharpList[_T] = self
        return ListEnumerator_1__ctor_3002E699(xs)

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None = None) -> IEnumerator[Any]:
        xs: FSharpList[_T] = self
        return get_enumerator(xs)


FSharpList_reflection = _expr68


def _expr69(gen0: TypeInfo) -> TypeInfo:
    return class_type("ListModule.ListEnumerator`1", [gen0], ListEnumerator_1)


class ListEnumerator_1(IDisposable, Generic[_T]):
    def __init__(self, xs: FSharpList[Any]) -> None:
        self.xs: FSharpList[_T] = xs
        self.it: FSharpList[_T] = self.xs
        self.current: _T = None

    def System_Collections_Generic_IEnumerator_1_get_Current(self, __unit: None = None) -> _T:
        _: ListEnumerator_1[_T] = self
        return _.current

    def System_Collections_IEnumerator_get_Current(self, __unit: None = None) -> Any:
        _: ListEnumerator_1[_T] = self
        return _.current

    def System_Collections_IEnumerator_MoveNext(self, __unit: None = None) -> bool:
        _: ListEnumerator_1[_T] = self
        match_value: FSharpList[_T] | None = _.it.tail
        if match_value is not None:
            t: FSharpList[_T] = match_value
            _.current = _.it.head
            _.it = t
            return True

        else:
            return False

    def System_Collections_IEnumerator_Reset(self, __unit: None = None) -> None:
        _: ListEnumerator_1[_T] = self
        _.it = _.xs
        _.current = None

    def Dispose(self, __unit: None = None) -> None:
        pass


ListEnumerator_1_reflection = _expr69


def ListEnumerator_1__ctor_3002E699(xs: FSharpList[Any]) -> ListEnumerator_1[_T]:
    return ListEnumerator_1(xs)


def FSharpList_get_Empty(__unit: None = None) -> FSharpList[Any]:
    return FSharpList(None, None)


def FSharpList_Cons_305B8EAC(x: _T, xs: FSharpList[_T]) -> FSharpList[_T]:
    return FSharpList(x, xs)


def FSharpList__get_IsEmpty(xs: FSharpList[Any]) -> bool:
    return xs.tail is None


def FSharpList__get_Length(xs: FSharpList[Any]) -> int:
    def loop(i_mut: int, xs_1_mut: FSharpList[_T], xs: Any = xs) -> int:
        while True:
            (i, xs_1) = (i_mut, xs_1_mut)
            match_value: FSharpList[_T] | None = xs_1.tail
            if match_value is not None:
                i_mut = i + 1
                xs_1_mut = match_value
                continue

            else:
                return i

            break

    return loop(0, xs)


def FSharpList__get_Head(xs: FSharpList[_T]) -> _T:
    match_value: FSharpList[_T] | None = xs.tail
    if match_value is not None:
        return xs.head

    else:
        raise Exception((SR_inputWasEmpty + "\\nParameter name: ") + "list")


def FSharpList__get_Tail(xs: FSharpList[_T]) -> FSharpList[_T]:
    match_value: FSharpList[_T] | None = xs.tail
    if match_value is not None:
        return match_value

    else:
        raise Exception((SR_inputWasEmpty + "\\nParameter name: ") + "list")


def FSharpList__get_Item_Z524259A4(xs: FSharpList[_T], index: int) -> _T:
    def loop(i_mut: int, xs_1_mut: FSharpList[_T], xs: Any = xs, index: Any = index) -> _T:
        while True:
            (i, xs_1) = (i_mut, xs_1_mut)
            match_value: FSharpList[_T] | None = xs_1.tail
            if match_value is not None:
                if i == index:
                    return xs_1.head

                else:
                    i_mut = i + 1
                    xs_1_mut = match_value
                    continue

            else:
                raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

            break

    return loop(0, xs)


def index_not_found(__unit: None = None) -> Any:
    raise Exception(SR_keyNotFoundAlt)


def empty(__unit: None = None) -> FSharpList[Any]:
    return FSharpList_get_Empty()


def cons(x: _T, xs: FSharpList[_T]) -> FSharpList[_T]:
    return FSharpList_Cons_305B8EAC(x, xs)


def singleton(x: __A | None = None) -> FSharpList[__A]:
    return FSharpList_Cons_305B8EAC(x, FSharpList_get_Empty())


def is_empty(xs: FSharpList[Any]) -> bool:
    return FSharpList__get_IsEmpty(xs)


def length(xs: FSharpList[Any]) -> int:
    return FSharpList__get_Length(xs)


def head(xs: FSharpList[_T]) -> _T:
    return FSharpList__get_Head(xs)


def try_head(xs: FSharpList[_T]) -> _T | None:
    if FSharpList__get_IsEmpty(xs):
        return None

    else:
        return some(FSharpList__get_Head(xs))


def tail(xs: FSharpList[_T]) -> FSharpList[_T]:
    return FSharpList__get_Tail(xs)


def try_last(xs_mut: FSharpList[_T]) -> _T | None:
    while True:
        (xs,) = (xs_mut,)
        if FSharpList__get_IsEmpty(xs):
            return None

        else:
            t: FSharpList[_T] = FSharpList__get_Tail(xs)
            if FSharpList__get_IsEmpty(t):
                return some(FSharpList__get_Head(xs))

            else:
                xs_mut = t
                continue

        break


def last(xs: FSharpList[_T]) -> _T:
    match_value: _T | None = try_last(xs)
    if match_value is None:
        raise Exception(SR_inputWasEmpty)

    else:
        return value_1(match_value)


def compare_with(comparer: Callable[[_T, _T], int], xs: FSharpList[_T], ys: FSharpList[_T]) -> int:
    def loop(
        xs_1_mut: FSharpList[_T], ys_1_mut: FSharpList[_T], comparer: Any = comparer, xs: Any = xs, ys: Any = ys
    ) -> int:
        while True:
            (xs_1, ys_1) = (xs_1_mut, ys_1_mut)
            matchValue: bool = FSharpList__get_IsEmpty(xs_1)
            matchValue_1: bool = FSharpList__get_IsEmpty(ys_1)
            if matchValue:
                if matchValue_1:
                    return 0

                else:
                    return -1

            elif matchValue_1:
                return 1

            else:
                c: int = comparer(FSharpList__get_Head(xs_1), FSharpList__get_Head(ys_1)) or 0
                if c == 0:
                    xs_1_mut = FSharpList__get_Tail(xs_1)
                    ys_1_mut = FSharpList__get_Tail(ys_1)
                    continue

                else:
                    return c

            break

    return loop(xs, ys)


def to_array(xs: FSharpList[_T]) -> Array[_T]:
    len_1: int = FSharpList__get_Length(xs) or 0
    res: Array[_T] = fill([0] * len_1, 0, len_1, None)

    def loop(i_mut: int, xs_1_mut: FSharpList[_T], xs: Any = xs) -> None:
        while True:
            (i, xs_1) = (i_mut, xs_1_mut)
            if not FSharpList__get_IsEmpty(xs_1):
                res[i] = FSharpList__get_Head(xs_1)
                i_mut = i + 1
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    loop(0, xs)
    return res


def fold(folder: Callable[[_STATE, _T], _STATE], state: _STATE, xs: FSharpList[_T]) -> _STATE:
    acc: _STATE = state
    xs_1: FSharpList[_T] = xs
    while not FSharpList__get_IsEmpty(xs_1):
        acc = folder(acc, head(xs_1))
        xs_1 = FSharpList__get_Tail(xs_1)
    return acc


def reverse(xs: FSharpList[_T]) -> FSharpList[_T]:
    def _arrow70(acc: FSharpList[_T], x: _T, xs: Any = xs) -> FSharpList[_T]:
        return FSharpList_Cons_305B8EAC(x, acc)

    return fold(_arrow70, FSharpList_get_Empty(), xs)


def fold_back(folder: Callable[[_T, _STATE], _STATE], xs: FSharpList[_T], state: _STATE) -> _STATE:
    return fold_back_1(folder, to_array(xs), state)


def fold_indexed(folder: Callable[[int, _STATE, _T], _STATE], state: _STATE, xs: FSharpList[_T]) -> _STATE:
    def loop(
        i_mut: int, acc_mut: _STATE, xs_1_mut: FSharpList[_T], folder: Any = folder, state: Any = state, xs: Any = xs
    ) -> _STATE:
        while True:
            (i, acc, xs_1) = (i_mut, acc_mut, xs_1_mut)
            if FSharpList__get_IsEmpty(xs_1):
                return acc

            else:
                i_mut = i + 1
                acc_mut = folder(i, acc, FSharpList__get_Head(xs_1))
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    return loop(0, state, xs)


def fold2(
    folder: Callable[[_STATE, _T1, _T2], _STATE], state: _STATE, xs: FSharpList[_T1], ys: FSharpList[_T2]
) -> _STATE:
    acc: _STATE = state
    xs_1: FSharpList[_T1] = xs
    ys_1: FSharpList[_T2] = ys
    while (not FSharpList__get_IsEmpty(ys_1)) if (not FSharpList__get_IsEmpty(xs_1)) else False:
        acc = folder(acc, FSharpList__get_Head(xs_1), FSharpList__get_Head(ys_1))
        xs_1 = FSharpList__get_Tail(xs_1)
        ys_1 = FSharpList__get_Tail(ys_1)
    return acc


def fold_back2(
    folder: Callable[[_T1, _T2, _STATE], _STATE], xs: FSharpList[_T1], ys: FSharpList[_T2], state: _STATE
) -> _STATE:
    return fold_back2_1(folder, to_array(xs), to_array(ys), state)


def unfold(gen: Callable[[_STATE], tuple[_T, _STATE] | None], state: _STATE) -> FSharpList[_T]:
    def loop(acc_mut: _STATE, node_mut: FSharpList[_T], gen: Any = gen, state: Any = state) -> FSharpList[_T]:
        while True:
            (acc, node) = (acc_mut, node_mut)
            match_value: tuple[_T, _STATE] | None = gen(acc)
            if match_value is not None:
                acc_mut = match_value[1]

                def _arrow71(__unit: None = None, acc: Any = acc, node: Any = node) -> FSharpList[_T]:
                    t: FSharpList[_T] = FSharpList(match_value[0], None)
                    node.tail = t
                    return t

                node_mut = _arrow71()
                continue

            else:
                return node

            break

    root: FSharpList[_T] = FSharpList_get_Empty()
    node_1: FSharpList[_T] = loop(state, root)
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    node_1.tail = t_2
    return FSharpList__get_Tail(root)


def iterate(action: Callable[[__A], None], xs: FSharpList[__A]) -> None:
    def _arrow72(unit_var: None, x: __A, action: Any = action, xs: Any = xs) -> None:
        action(x)

    fold(_arrow72, None, xs)


def iterate2(action: Callable[[__A, __B], None], xs: FSharpList[__A], ys: FSharpList[__B]) -> None:
    def _arrow73(unit_var: None, x: __A, y: __B, action: Any = action, xs: Any = xs, ys: Any = ys) -> None:
        action(x, y)

    fold2(_arrow73, None, xs, ys)


def iterate_indexed(action: Callable[[int, __A], None], xs: FSharpList[__A]) -> None:
    def _arrow74(i: int, x: __A, action: Any = action, xs: Any = xs) -> int:
        action(i, x)
        return i + 1

    ignore(fold(_arrow74, 0, xs))


def iterate_indexed2(action: Callable[[int, __A, __B], None], xs: FSharpList[__A], ys: FSharpList[__B]) -> None:
    def _arrow75(i: int, x: __A, y: __B, action: Any = action, xs: Any = xs, ys: Any = ys) -> int:
        action(i, x, y)
        return i + 1

    ignore(fold2(_arrow75, 0, xs, ys))


def to_seq(xs: FSharpList[_T]) -> IEnumerable_1[_T]:
    return xs


def of_array_with_tail(xs: Array[_T], tail_1: FSharpList[_T]) -> FSharpList[_T]:
    res: FSharpList[_T] = tail_1
    for i in range(len(xs) - 1, 0 - 1, -1):
        res = FSharpList_Cons_305B8EAC(xs[i], res)
    return res


def of_array(xs: Array[_T]) -> FSharpList[_T]:
    return of_array_with_tail(xs, FSharpList_get_Empty())


def of_seq(xs: IEnumerable_1[_T]) -> FSharpList[_T]:
    if is_array_like(xs):
        return of_array(xs)

    elif isinstance(xs, FSharpList):
        return xs

    else:
        root: FSharpList[_T] = FSharpList_get_Empty()
        node: FSharpList[_T] = root
        with get_enumerator(xs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                x: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()

                def _arrow76(__unit: None = None, xs: Any = xs) -> FSharpList[_T]:
                    xs_3: FSharpList[_T] = node
                    t: FSharpList[_T] = FSharpList(x, None)
                    xs_3.tail = t
                    return t

                node = _arrow76()
        xs_5: FSharpList[_T] = node
        t_2: FSharpList[_T] = FSharpList_get_Empty()
        xs_5.tail = t_2
        return FSharpList__get_Tail(root)


def concat(lists: IEnumerable_1[FSharpList[_T]]) -> FSharpList[_T]:
    root: FSharpList[_T] = FSharpList_get_Empty()
    node: FSharpList[_T] = root

    def action(xs: FSharpList[_T], lists: Any = lists) -> None:
        nonlocal node

        def _arrow77(acc: FSharpList[_T], x: _T, xs: Any = xs) -> FSharpList[_T]:
            t: FSharpList[_T] = FSharpList(x, None)
            acc.tail = t
            return t

        node = fold(_arrow77, node, xs)

    if is_array_like(lists):
        iterate_1(action, lists)

    elif isinstance(lists, FSharpList):
        iterate(action, lists)

    else:
        with get_enumerator(lists) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                action(enumerator.System_Collections_Generic_IEnumerator_1_get_Current())

    xs_6: FSharpList[_T] = node
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    xs_6.tail = t_2
    return FSharpList__get_Tail(root)


def scan(folder: Callable[[_STATE, _T], _STATE], state: _STATE, xs: FSharpList[_T]) -> FSharpList[_STATE]:
    root: FSharpList[_STATE] = FSharpList_get_Empty()
    node: FSharpList[_STATE]
    t: FSharpList[_STATE] = FSharpList(state, None)
    root.tail = t
    node = t
    acc: _STATE = state
    xs_3: FSharpList[_T] = xs
    while not FSharpList__get_IsEmpty(xs_3):
        acc = folder(acc, FSharpList__get_Head(xs_3))

        def _arrow78(__unit: None = None, folder: Any = folder, state: Any = state, xs: Any = xs) -> FSharpList[_STATE]:
            xs_4: FSharpList[_STATE] = node
            t_2: FSharpList[_STATE] = FSharpList(acc, None)
            xs_4.tail = t_2
            return t_2

        node = _arrow78()
        xs_3 = FSharpList__get_Tail(xs_3)
    xs_6: FSharpList[_STATE] = node
    t_4: FSharpList[_STATE] = FSharpList_get_Empty()
    xs_6.tail = t_4
    return FSharpList__get_Tail(root)


def scan_back(folder: Callable[[_T, _STATE], _STATE], xs: FSharpList[_T], state: _STATE) -> FSharpList[_STATE]:
    return of_array(scan_back_1(folder, to_array(xs), state, None))


def append(xs: FSharpList[_T], ys: FSharpList[_T]) -> FSharpList[_T]:
    def _arrow79(acc: FSharpList[_T], x: _T, xs: Any = xs, ys: Any = ys) -> FSharpList[_T]:
        return FSharpList_Cons_305B8EAC(x, acc)

    return fold(_arrow79, ys, reverse(xs))


def collect(mapping: Callable[[_T], FSharpList[_U]], xs: FSharpList[_T]) -> FSharpList[_U]:
    root: FSharpList[_U] = FSharpList_get_Empty()
    node: FSharpList[_U] = root
    ys: FSharpList[_T] = xs
    while not FSharpList__get_IsEmpty(ys):
        zs: FSharpList[_U] = mapping(FSharpList__get_Head(ys))
        while not FSharpList__get_IsEmpty(zs):

            def _arrow80(__unit: None = None, mapping: Any = mapping, xs: Any = xs) -> FSharpList[_U]:
                xs_1: FSharpList[_U] = node
                t: FSharpList[_U] = FSharpList(FSharpList__get_Head(zs), None)
                xs_1.tail = t
                return t

            node = _arrow80()
            zs = FSharpList__get_Tail(zs)
        ys = FSharpList__get_Tail(ys)
    xs_3: FSharpList[_U] = node
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    xs_3.tail = t_2
    return FSharpList__get_Tail(root)


def map_indexed(mapping: Callable[[int, _T], _U], xs: FSharpList[_T]) -> FSharpList[_U]:
    root: FSharpList[_U] = FSharpList_get_Empty()

    def folder(i: int, acc: FSharpList[_U], x: _T, mapping: Any = mapping, xs: Any = xs) -> FSharpList[_U]:
        t: FSharpList[_U] = FSharpList(mapping(i, x), None)
        acc.tail = t
        return t

    node: FSharpList[_U] = fold_indexed(folder, root, xs)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def map(mapping: Callable[[_T], _U], xs: FSharpList[_T]) -> FSharpList[_U]:
    root: FSharpList[_U] = FSharpList_get_Empty()

    def folder(acc: FSharpList[_U], x: _T, mapping: Any = mapping, xs: Any = xs) -> FSharpList[_U]:
        t: FSharpList[_U] = FSharpList(mapping(x), None)
        acc.tail = t
        return t

    node: FSharpList[_U] = fold(folder, root, xs)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def indexed(xs: FSharpList[__A]) -> FSharpList[tuple[int, __A]]:
    def _arrow81(i: int, x: __A, xs: Any = xs) -> tuple[int, __A]:
        return (i, x)

    return map_indexed(_arrow81, xs)


def map2(mapping: Callable[[_T1, _T2], _U], xs: FSharpList[_T1], ys: FSharpList[_T2]) -> FSharpList[_U]:
    root: FSharpList[_U] = FSharpList_get_Empty()

    def folder(
        acc: FSharpList[_U], x: _T1, y: _T2, mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> FSharpList[_U]:
        t: FSharpList[_U] = FSharpList(mapping(x, y), None)
        acc.tail = t
        return t

    node: FSharpList[_U] = fold2(folder, root, xs, ys)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def map_indexed2(mapping: Callable[[int, _T1, _T2], _U], xs: FSharpList[_T1], ys: FSharpList[_T2]) -> FSharpList[_U]:
    def loop(
        i_mut: int,
        acc_mut: FSharpList[_U],
        xs_1_mut: FSharpList[_T1],
        ys_1_mut: FSharpList[_T2],
        mapping: Any = mapping,
        xs: Any = xs,
        ys: Any = ys,
    ) -> FSharpList[_U]:
        while True:
            (i, acc, xs_1, ys_1) = (i_mut, acc_mut, xs_1_mut, ys_1_mut)
            if True if FSharpList__get_IsEmpty(xs_1) else FSharpList__get_IsEmpty(ys_1):
                return acc

            else:
                i_mut = i + 1

                def _arrow82(
                    __unit: None = None, i: Any = i, acc: Any = acc, xs_1: Any = xs_1, ys_1: Any = ys_1
                ) -> FSharpList[_U]:
                    t: FSharpList[_U] = FSharpList(
                        mapping(i, FSharpList__get_Head(xs_1), FSharpList__get_Head(ys_1)), None
                    )
                    acc.tail = t
                    return t

                acc_mut = _arrow82()
                xs_1_mut = FSharpList__get_Tail(xs_1)
                ys_1_mut = FSharpList__get_Tail(ys_1)
                continue

            break

    root: FSharpList[_U] = FSharpList_get_Empty()
    node_1: FSharpList[_U] = loop(0, root, xs, ys)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node_1.tail = t_2
    return FSharpList__get_Tail(root)


def map3(
    mapping: Callable[[_T1, _T2, _T3], _U], xs: FSharpList[_T1], ys: FSharpList[_T2], zs: FSharpList[_T3]
) -> FSharpList[_U]:
    def loop(
        acc_mut: FSharpList[_U],
        xs_1_mut: FSharpList[_T1],
        ys_1_mut: FSharpList[_T2],
        zs_1_mut: FSharpList[_T3],
        mapping: Any = mapping,
        xs: Any = xs,
        ys: Any = ys,
        zs: Any = zs,
    ) -> FSharpList[_U]:
        while True:
            (acc, xs_1, ys_1, zs_1) = (acc_mut, xs_1_mut, ys_1_mut, zs_1_mut)
            if (
                True
                if (True if FSharpList__get_IsEmpty(xs_1) else FSharpList__get_IsEmpty(ys_1))
                else FSharpList__get_IsEmpty(zs_1)
            ):
                return acc

            else:

                def _arrow83(
                    __unit: None = None, acc: Any = acc, xs_1: Any = xs_1, ys_1: Any = ys_1, zs_1: Any = zs_1
                ) -> FSharpList[_U]:
                    t: FSharpList[_U] = FSharpList(
                        mapping(FSharpList__get_Head(xs_1), FSharpList__get_Head(ys_1), FSharpList__get_Head(zs_1)),
                        None,
                    )
                    acc.tail = t
                    return t

                acc_mut = _arrow83()
                xs_1_mut = FSharpList__get_Tail(xs_1)
                ys_1_mut = FSharpList__get_Tail(ys_1)
                zs_1_mut = FSharpList__get_Tail(zs_1)
                continue

            break

    root: FSharpList[_U] = FSharpList_get_Empty()
    node_1: FSharpList[_U] = loop(root, xs, ys, zs)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node_1.tail = t_2
    return FSharpList__get_Tail(root)


def map_fold(
    mapping: Callable[[_STATE, _T], tuple[_RESULT, _STATE]], state: _STATE, xs: FSharpList[_T]
) -> tuple[FSharpList[_RESULT], _STATE]:
    root: FSharpList[_RESULT] = FSharpList_get_Empty()

    def folder(
        tupled_arg: tuple[FSharpList[_RESULT], _STATE], x: _T, mapping: Any = mapping, state: Any = state, xs: Any = xs
    ) -> tuple[FSharpList[_RESULT], _STATE]:
        pattern_input: tuple[_RESULT, _STATE] = mapping(tupled_arg[1], x)

        def _arrow84(__unit: None = None, tupled_arg: Any = tupled_arg, x: Any = x) -> FSharpList[_RESULT]:
            t: FSharpList[_RESULT] = FSharpList(pattern_input[0], None)
            tupled_arg[0].tail = t
            return t

        return (_arrow84(), pattern_input[1])

    pattern_input_1: tuple[FSharpList[_RESULT], _STATE] = fold(folder, (root, state), xs)
    t_2: FSharpList[_RESULT] = FSharpList_get_Empty()
    pattern_input_1[0].tail = t_2
    return (FSharpList__get_Tail(root), pattern_input_1[1])


def map_fold_back(
    mapping: Callable[[_T, _STATE], tuple[_RESULT, _STATE]], xs: FSharpList[_T], state: _STATE
) -> tuple[FSharpList[_RESULT], _STATE]:
    def _arrow85(
        acc: _STATE, x: _T, mapping: Any = mapping, xs: Any = xs, state: Any = state
    ) -> tuple[_RESULT, _STATE]:
        return mapping(x, acc)

    return map_fold(_arrow85, state, reverse(xs))


def try_pick(f: Callable[[_T], __A | None], xs: FSharpList[_T]) -> __A | None:
    def loop(xs_1_mut: FSharpList[_T], f: Any = f, xs: Any = xs) -> __A | None:
        while True:
            (xs_1,) = (xs_1_mut,)
            if FSharpList__get_IsEmpty(xs_1):
                return None

            else:
                match_value: __A | None = f(FSharpList__get_Head(xs_1))
                if match_value is None:
                    xs_1_mut = FSharpList__get_Tail(xs_1)
                    continue

                else:
                    return match_value

            break

    return loop(xs)


def pick(f: Callable[[__A], __B | None], xs: FSharpList[__A]) -> __B:
    match_value: __B | None = try_pick(f, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def try_find(f: Callable[[__A], bool], xs: FSharpList[__A]) -> __A | None:
    def _arrow86(x: __A | None = None, f: Any = f, xs: Any = xs) -> __A | None:
        return some(x) if f(x) else None

    return try_pick(_arrow86, xs)


def find(f: Callable[[__A], bool], xs: FSharpList[__A]) -> __A:
    match_value: __A | None = try_find(f, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def try_find_back(f: Callable[[__A], bool], xs: FSharpList[__A]) -> __A | None:
    return try_find_back_1(f, to_array(xs))


def find_back(f: Callable[[__A], bool], xs: FSharpList[__A]) -> __A:
    match_value: __A | None = try_find_back(f, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def try_find_index(f: Callable[[_T], bool], xs: FSharpList[_T]) -> int | None:
    def loop(i_mut: int, xs_1_mut: FSharpList[_T], f: Any = f, xs: Any = xs) -> int | None:
        while True:
            (i, xs_1) = (i_mut, xs_1_mut)
            if FSharpList__get_IsEmpty(xs_1):
                return None

            elif f(FSharpList__get_Head(xs_1)):
                return i

            else:
                i_mut = i + 1
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    return loop(0, xs)


def find_index(f: Callable[[__A], bool], xs: FSharpList[__A]) -> int:
    match_value: int | None = try_find_index(f, xs)
    if match_value is None:
        index_not_found()
        return -1

    else:
        return match_value


def try_find_index_back(f: Callable[[__A], bool], xs: FSharpList[__A]) -> int | None:
    return try_find_index_back_1(f, to_array(xs))


def find_index_back(f: Callable[[__A], bool], xs: FSharpList[__A]) -> int:
    match_value: int | None = try_find_index_back(f, xs)
    if match_value is None:
        index_not_found()
        return -1

    else:
        return match_value


def try_item(n: int, xs: FSharpList[_T]) -> _T | None:
    def loop(i_mut: int, xs_1_mut: FSharpList[_T], n: Any = n, xs: Any = xs) -> _T | None:
        while True:
            (i, xs_1) = (i_mut, xs_1_mut)
            if FSharpList__get_IsEmpty(xs_1):
                return None

            elif i == n:
                return some(FSharpList__get_Head(xs_1))

            else:
                i_mut = i + 1
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    return loop(0, xs)


def item(n: int, xs: FSharpList[_T]) -> _T:
    return FSharpList__get_Item_Z524259A4(xs, n)


def filter(f: Callable[[_T], bool], xs: FSharpList[_T]) -> FSharpList[_T]:
    root: FSharpList[_T] = FSharpList_get_Empty()

    def folder(acc: FSharpList[_T], x: _T, f: Any = f, xs: Any = xs) -> FSharpList[_T]:
        if f(x):
            t: FSharpList[_T] = FSharpList(x, None)
            acc.tail = t
            return t

        else:
            return acc

    node: FSharpList[_T] = fold(folder, root, xs)
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def partition(f: Callable[[_T], bool], xs: FSharpList[_T]) -> tuple[FSharpList[_T], FSharpList[_T]]:
    matchValue: FSharpList[_T] = FSharpList_get_Empty()
    root2: FSharpList[_T] = FSharpList_get_Empty()
    root1: FSharpList[_T] = matchValue

    def folder(
        tupled_arg: tuple[FSharpList[_T], FSharpList[_T]], x: _T, f: Any = f, xs: Any = xs
    ) -> tuple[FSharpList[_T], FSharpList[_T]]:
        lacc: FSharpList[_T] = tupled_arg[0]
        racc: FSharpList[_T] = tupled_arg[1]
        if f(x):

            def _arrow87(__unit: None = None, tupled_arg: Any = tupled_arg, x: Any = x) -> FSharpList[_T]:
                t: FSharpList[_T] = FSharpList(x, None)
                lacc.tail = t
                return t

            return (_arrow87(), racc)

        else:

            def _arrow88(__unit: None = None, tupled_arg: Any = tupled_arg, x: Any = x) -> FSharpList[_T]:
                t_2: FSharpList[_T] = FSharpList(x, None)
                racc.tail = t_2
                return t_2

            return (lacc, _arrow88())

    pattern_input_1: tuple[FSharpList[_T], FSharpList[_T]] = fold(folder, (root1, root2), xs)
    t_4: FSharpList[_T] = FSharpList_get_Empty()
    pattern_input_1[0].tail = t_4
    t_5: FSharpList[_T] = FSharpList_get_Empty()
    pattern_input_1[1].tail = t_5
    return (FSharpList__get_Tail(root1), FSharpList__get_Tail(root2))


def choose(f: Callable[[_T], _U | None], xs: FSharpList[_T]) -> FSharpList[_U]:
    root: FSharpList[_U] = FSharpList_get_Empty()

    def folder(acc: FSharpList[_U], x: _T, f: Any = f, xs: Any = xs) -> FSharpList[_U]:
        match_value: _U | None = f(x)
        if match_value is None:
            return acc

        else:
            t: FSharpList[_U] = FSharpList(value_1(match_value), None)
            acc.tail = t
            return t

    node: FSharpList[_U] = fold(folder, root, xs)
    t_2: FSharpList[_U] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def contains(value: _T, xs: FSharpList[_T], eq: IEqualityComparer_1[Any]) -> bool:
    def _arrow89(v: _T | None = None, value: Any = value, xs: Any = xs, eq: Any = eq) -> bool:
        return eq.Equals(value, v)

    return try_find_index(_arrow89, xs) is not None


def initialize(n: int, f: Callable[[int], _T]) -> FSharpList[_T]:
    root: FSharpList[_T] = FSharpList_get_Empty()
    node: FSharpList[_T] = root
    for i in range(0, (n - 1) + 1, 1):

        def _arrow90(__unit: None = None, n: Any = n, f: Any = f) -> FSharpList[_T]:
            xs: FSharpList[_T] = node
            t: FSharpList[_T] = FSharpList(f(i), None)
            xs.tail = t
            return t

        node = _arrow90()
    xs_2: FSharpList[_T] = node
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    xs_2.tail = t_2
    return FSharpList__get_Tail(root)


def replicate(n: int, x: __A) -> FSharpList[__A]:
    def _arrow91(_arg: int, n: Any = n, x: Any = x) -> __A:
        return x

    return initialize(n, _arrow91)


def reduce(f: Callable[[_T, _T], _T], xs: FSharpList[_T]) -> _T:
    if FSharpList__get_IsEmpty(xs):
        raise Exception(SR_inputWasEmpty)

    else:
        return fold(f, head(xs), tail(xs))


def reduce_back(f: Callable[[_T, _T], _T], xs: FSharpList[_T]) -> _T:
    if FSharpList__get_IsEmpty(xs):
        raise Exception(SR_inputWasEmpty)

    else:
        return fold_back(f, tail(xs), head(xs))


def for_all(f: Callable[[__A], bool], xs: FSharpList[__A]) -> bool:
    def _arrow92(acc: bool, x: __A, f: Any = f, xs: Any = xs) -> bool:
        return f(x) if acc else False

    return fold(_arrow92, True, xs)


def for_all2(f: Callable[[__A, __B], bool], xs: FSharpList[__A], ys: FSharpList[__B]) -> bool:
    def _arrow93(acc: bool, x: __A, y: __B, f: Any = f, xs: Any = xs, ys: Any = ys) -> bool:
        return f(x, y) if acc else False

    return fold2(_arrow93, True, xs, ys)


def exists(f: Callable[[__A], bool], xs: FSharpList[__A]) -> bool:
    return try_find_index(f, xs) is not None


def exists2(f_mut: Callable[[_T1, _T2], bool], xs_mut: FSharpList[_T1], ys_mut: FSharpList[_T2]) -> bool:
    while True:
        (f, xs, ys) = (f_mut, xs_mut, ys_mut)
        matchValue: bool = FSharpList__get_IsEmpty(xs)
        matchValue_1: bool = FSharpList__get_IsEmpty(ys)
        (pattern_matching_result,) = (None,)
        if matchValue:
            if matchValue_1:
                pattern_matching_result = 0

            else:
                pattern_matching_result = 2

        elif matchValue_1:
            pattern_matching_result = 2

        else:
            pattern_matching_result = 1

        if pattern_matching_result == 0:
            return False

        elif pattern_matching_result == 1:
            if f(FSharpList__get_Head(xs), FSharpList__get_Head(ys)):
                return True

            else:
                f_mut = f
                xs_mut = FSharpList__get_Tail(xs)
                ys_mut = FSharpList__get_Tail(ys)
                continue

        elif pattern_matching_result == 2:
            raise Exception((SR_differentLengths + "\\nParameter name: ") + "list2")

        break


def unzip(xs: FSharpList[tuple[__A, __B]]) -> tuple[FSharpList[__A], FSharpList[__B]]:
    def _arrow94(
        tupled_arg: tuple[__A, __B], tupled_arg_1: tuple[FSharpList[__A], FSharpList[__B]], xs: Any = xs
    ) -> tuple[FSharpList[__A], FSharpList[__B]]:
        return (
            FSharpList_Cons_305B8EAC(tupled_arg[0], tupled_arg_1[0]),
            FSharpList_Cons_305B8EAC(tupled_arg[1], tupled_arg_1[1]),
        )

    return fold_back(_arrow94, xs, (FSharpList_get_Empty(), FSharpList_get_Empty()))


def unzip3(xs: FSharpList[tuple[__A, __B, __C]]) -> tuple[FSharpList[__A], FSharpList[__B], FSharpList[__C]]:
    def _arrow95(
        tupled_arg: tuple[__A, __B, __C],
        tupled_arg_1: tuple[FSharpList[__A], FSharpList[__B], FSharpList[__C]],
        xs: Any = xs,
    ) -> tuple[FSharpList[__A], FSharpList[__B], FSharpList[__C]]:
        return (
            FSharpList_Cons_305B8EAC(tupled_arg[0], tupled_arg_1[0]),
            FSharpList_Cons_305B8EAC(tupled_arg[1], tupled_arg_1[1]),
            FSharpList_Cons_305B8EAC(tupled_arg[2], tupled_arg_1[2]),
        )

    return fold_back(_arrow95, xs, (FSharpList_get_Empty(), FSharpList_get_Empty(), FSharpList_get_Empty()))


def zip(xs: FSharpList[__A], ys: FSharpList[__B]) -> FSharpList[tuple[__A, __B]]:
    def _arrow96(x: __A, y: __B, xs: Any = xs, ys: Any = ys) -> tuple[__A, __B]:
        return (x, y)

    return map2(_arrow96, xs, ys)


def zip3(xs: FSharpList[__A], ys: FSharpList[__B], zs: FSharpList[__C]) -> FSharpList[tuple[__A, __B, __C]]:
    def _arrow97(x: __A, y: __B, z: __C, xs: Any = xs, ys: Any = ys, zs: Any = zs) -> tuple[__A, __B, __C]:
        return (x, y, z)

    return map3(_arrow97, xs, ys, zs)


def sort_with(comparer: Callable[[_T, _T], int], xs: FSharpList[_T]) -> FSharpList[_T]:
    arr: Array[_T] = to_array(xs)
    sort_in_place_with(comparer, arr)
    return of_array(arr)


def sort(xs: FSharpList[_T], comparer: IComparer_1[_T]) -> FSharpList[_T]:
    def _arrow98(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(x, y)

    return sort_with(_arrow98, xs)


def sort_by(projection: Callable[[_T], _U], xs: FSharpList[_T], comparer: IComparer_1[_U]) -> FSharpList[_T]:
    def _arrow99(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(projection(x), projection(y))

    return sort_with(_arrow99, xs)


def sort_descending(xs: FSharpList[_T], comparer: IComparer_1[_T]) -> FSharpList[_T]:
    def _arrow100(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(x, y) * -1

    return sort_with(_arrow100, xs)


def sort_by_descending(projection: Callable[[_T], _U], xs: FSharpList[_T], comparer: IComparer_1[_U]) -> FSharpList[_T]:
    def _arrow101(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(projection(x), projection(y)) * -1

    return sort_with(_arrow101, xs)


def sum(xs: FSharpList[_T], adder: IGenericAdder_1[_T]) -> _T:
    def _arrow102(acc: _T, x: _T, xs: Any = xs, adder: Any = adder) -> _T:
        return adder.Add(acc, x)

    return fold(_arrow102, adder.GetZero(), xs)


def sum_by(f: Callable[[_T], _U], xs: FSharpList[_T], adder: IGenericAdder_1[_U]) -> _U:
    def _arrow103(acc: _U, x: _T, f: Any = f, xs: Any = xs, adder: Any = adder) -> _U:
        return adder.Add(acc, f(x))

    return fold(_arrow103, adder.GetZero(), xs)


def max_by(projection: Callable[[_T], _U], xs: FSharpList[_T], comparer: IComparer_1[_U]) -> _T:
    def _arrow104(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> _T:
        return y if (comparer.Compare(projection(y), projection(x)) > 0) else x

    return reduce(_arrow104, xs)


def max(xs: FSharpList[_T], comparer: IComparer_1[_T]) -> _T:
    def _arrow105(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> _T:
        return y if (comparer.Compare(y, x) > 0) else x

    return reduce(_arrow105, xs)


def min_by(projection: Callable[[_T], _U], xs: FSharpList[_T], comparer: IComparer_1[_U]) -> _T:
    def _arrow106(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> _T:
        return x if (comparer.Compare(projection(y), projection(x)) > 0) else y

    return reduce(_arrow106, xs)


def min(xs: FSharpList[_T], comparer: IComparer_1[_T]) -> _T:
    def _arrow107(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> _T:
        return x if (comparer.Compare(y, x) > 0) else y

    return reduce(_arrow107, xs)


def average(xs: FSharpList[_T], averager: IGenericAverager_1[_T]) -> _T:
    count: int = 0

    def folder(acc: _T, x: _T, xs: Any = xs, averager: Any = averager) -> _T:
        nonlocal count
        count = (count + 1) or 0
        return averager.Add(acc, x)

    total: _T = fold(folder, averager.GetZero(), xs)
    return averager.DivideByInt(total, count)


def average_by(f: Callable[[_T], _U], xs: FSharpList[_T], averager: IGenericAverager_1[_U]) -> _U:
    count: int = 0

    def _arrow108(acc: _U, x: _T, f: Any = f, xs: Any = xs, averager: Any = averager) -> _U:
        nonlocal count
        count = (count + 1) or 0
        return averager.Add(acc, f(x))

    total: _U = fold(_arrow108, averager.GetZero(), xs)
    return averager.DivideByInt(total, count)


def permute(f: Callable[[int], int], xs: FSharpList[_T]) -> FSharpList[_T]:
    return of_array(permute_1(f, to_array(xs)))


def chunk_by_size(chunk_size: int, xs: FSharpList[_T]) -> FSharpList[FSharpList[_T]]:
    def mapping(xs_1: Array[_T], chunk_size: Any = chunk_size, xs: Any = xs) -> FSharpList[_T]:
        return of_array(xs_1)

    return of_array(map_1(mapping, chunk_by_size_1(chunk_size, to_array(xs)), None))


def all_pairs(xs: FSharpList[_T1], ys: FSharpList[_T2]) -> FSharpList[tuple[_T1, _T2]]:
    root: FSharpList[tuple[_T1, _T2]] = FSharpList_get_Empty()
    node: FSharpList[tuple[_T1, _T2]] = root

    def _arrow111(x: _T1 | None = None, xs: Any = xs, ys: Any = ys) -> None:
        def _arrow110(y: _T2 | None = None) -> None:
            nonlocal node

            def _arrow109(__unit: None = None) -> FSharpList[tuple[_T1, _T2]]:
                xs_1: FSharpList[tuple[_T1, _T2]] = node
                t: FSharpList[tuple[_T1, _T2]] = FSharpList((x, y), None)
                xs_1.tail = t
                return t

            node = _arrow109()

        iterate(_arrow110, ys)

    iterate(_arrow111, xs)
    xs_3: FSharpList[tuple[_T1, _T2]] = node
    t_2: FSharpList[tuple[_T1, _T2]] = FSharpList_get_Empty()
    xs_3.tail = t_2
    return FSharpList__get_Tail(root)


def skip(count_mut: int, xs_mut: FSharpList[_T]) -> FSharpList[_T]:
    while True:
        (count, xs) = (count_mut, xs_mut)
        if count <= 0:
            return xs

        elif FSharpList__get_IsEmpty(xs):
            raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "list")

        else:
            count_mut = count - 1
            xs_mut = FSharpList__get_Tail(xs)
            continue

        break


def skip_while(predicate_mut: Callable[[_T], bool], xs_mut: FSharpList[_T]) -> FSharpList[_T]:
    while True:
        (predicate, xs) = (predicate_mut, xs_mut)
        if FSharpList__get_IsEmpty(xs):
            return xs

        elif not predicate(FSharpList__get_Head(xs)):
            return xs

        else:
            predicate_mut = predicate
            xs_mut = FSharpList__get_Tail(xs)
            continue

        break


def take(count: int, xs: FSharpList[_T]) -> FSharpList[_T]:
    if count < 0:
        raise Exception((SR_inputMustBeNonNegative + "\\nParameter name: ") + "count")

    def loop(
        i_mut: int, acc_mut: FSharpList[_T], xs_1_mut: FSharpList[_T], count: Any = count, xs: Any = xs
    ) -> FSharpList[_T]:
        while True:
            (i, acc, xs_1) = (i_mut, acc_mut, xs_1_mut)
            if i <= 0:
                return acc

            elif FSharpList__get_IsEmpty(xs_1):
                raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "list")

            else:
                i_mut = i - 1

                def _arrow112(__unit: None = None, i: Any = i, acc: Any = acc, xs_1: Any = xs_1) -> FSharpList[_T]:
                    t: FSharpList[_T] = FSharpList(FSharpList__get_Head(xs_1), None)
                    acc.tail = t
                    return t

                acc_mut = _arrow112()
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    root: FSharpList[_T] = FSharpList_get_Empty()
    node: FSharpList[_T] = loop(count, root, xs)
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def take_while(predicate: Callable[[_T], bool], xs: FSharpList[_T]) -> FSharpList[_T]:
    def loop(
        acc_mut: FSharpList[_T], xs_1_mut: FSharpList[_T], predicate: Any = predicate, xs: Any = xs
    ) -> FSharpList[_T]:
        while True:
            (acc, xs_1) = (acc_mut, xs_1_mut)
            if FSharpList__get_IsEmpty(xs_1):
                return acc

            elif not predicate(FSharpList__get_Head(xs_1)):
                return acc

            else:

                def _arrow113(__unit: None = None, acc: Any = acc, xs_1: Any = xs_1) -> FSharpList[_T]:
                    t: FSharpList[_T] = FSharpList(FSharpList__get_Head(xs_1), None)
                    acc.tail = t
                    return t

                acc_mut = _arrow113()
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    root: FSharpList[_T] = FSharpList_get_Empty()
    node: FSharpList[_T] = loop(root, xs)
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def truncate(count: int, xs: FSharpList[_T]) -> FSharpList[_T]:
    def loop(
        i_mut: int, acc_mut: FSharpList[_T], xs_1_mut: FSharpList[_T], count: Any = count, xs: Any = xs
    ) -> FSharpList[_T]:
        while True:
            (i, acc, xs_1) = (i_mut, acc_mut, xs_1_mut)
            if i <= 0:
                return acc

            elif FSharpList__get_IsEmpty(xs_1):
                return acc

            else:
                i_mut = i - 1

                def _arrow114(__unit: None = None, i: Any = i, acc: Any = acc, xs_1: Any = xs_1) -> FSharpList[_T]:
                    t: FSharpList[_T] = FSharpList(FSharpList__get_Head(xs_1), None)
                    acc.tail = t
                    return t

                acc_mut = _arrow114()
                xs_1_mut = FSharpList__get_Tail(xs_1)
                continue

            break

    root: FSharpList[_T] = FSharpList_get_Empty()
    node: FSharpList[_T] = loop(count, root, xs)
    t_2: FSharpList[_T] = FSharpList_get_Empty()
    node.tail = t_2
    return FSharpList__get_Tail(root)


def get_slice(start_index: int | None, end_index: int | None, xs: FSharpList[_T]) -> FSharpList[_T]:
    len_1: int = length(xs) or 0
    start_index_1: int
    index: int = default_arg(start_index, 0) or 0
    start_index_1 = 0 if (index < 0) else index
    end_index_1: int
    index_1: int = default_arg(end_index, len_1 - 1) or 0
    end_index_1 = (len_1 - 1) if (index_1 >= len_1) else index_1
    if end_index_1 < start_index_1:
        return FSharpList_get_Empty()

    else:
        return take((end_index_1 - start_index_1) + 1, skip(start_index_1, xs))


def split_at(index: int, xs: FSharpList[_T]) -> tuple[FSharpList[_T], FSharpList[_T]]:
    if index < 0:
        raise Exception((SR_inputMustBeNonNegative + "\\nParameter name: ") + "index")

    if index > FSharpList__get_Length(xs):
        raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "index")

    return (take(index, xs), skip(index, xs))


def exactly_one(xs: FSharpList[_T]) -> _T:
    if FSharpList__get_IsEmpty(xs):
        raise Exception((SR_inputSequenceEmpty + "\\nParameter name: ") + "list")

    elif FSharpList__get_IsEmpty(FSharpList__get_Tail(xs)):
        return FSharpList__get_Head(xs)

    else:
        raise Exception((SR_inputSequenceTooLong + "\\nParameter name: ") + "list")


def try_exactly_one(xs: FSharpList[_T]) -> _T | None:
    if FSharpList__get_IsEmpty(FSharpList__get_Tail(xs)) if (not FSharpList__get_IsEmpty(xs)) else False:
        return some(FSharpList__get_Head(xs))

    else:
        return None


def where(predicate: Callable[[_T], bool], xs: FSharpList[_T]) -> FSharpList[_T]:
    return filter(predicate, xs)


def pairwise(xs: FSharpList[_T]) -> FSharpList[tuple[_T, _T]]:
    return of_array(pairwise_1(to_array(xs)))


def windowed(window_size: int, xs: FSharpList[_T]) -> FSharpList[FSharpList[_T]]:
    def mapping(xs_1: Array[_T], window_size: Any = window_size, xs: Any = xs) -> FSharpList[_T]:
        return of_array(xs_1)

    return of_array(map_1(mapping, windowed_1(window_size, to_array(xs)), None))


def split_into(chunks: int, xs: FSharpList[_T]) -> FSharpList[FSharpList[_T]]:
    def mapping(xs_1: Array[_T], chunks: Any = chunks, xs: Any = xs) -> FSharpList[_T]:
        return of_array(xs_1)

    return of_array(map_1(mapping, split_into_1(chunks, to_array(xs)), None))


def transpose(lists: IEnumerable_1[FSharpList[_T]]) -> FSharpList[FSharpList[_T]]:
    def mapping_1(xs_1: Array[_T], lists: Any = lists) -> FSharpList[_T]:
        return of_array(xs_1)

    def mapping(xs: FSharpList[_T], lists: Any = lists) -> Array[_T]:
        return to_array(xs)

    return of_array(map_1(mapping_1, transpose_1(map_1(mapping, list(lists), None), None), None))


def insert_at(index: int, y: _T, xs: FSharpList[_T]) -> FSharpList[_T]:
    i: int = -1
    is_done: bool = False

    def folder(acc: FSharpList[_T], x: _T, index: Any = index, y: Any = y, xs: Any = xs) -> FSharpList[_T]:
        nonlocal i, is_done
        i = (i + 1) or 0
        if i == index:
            is_done = True
            return FSharpList_Cons_305B8EAC(x, FSharpList_Cons_305B8EAC(y, acc))

        else:
            return FSharpList_Cons_305B8EAC(x, acc)

    result: FSharpList[_T] = fold(folder, FSharpList_get_Empty(), xs)

    def _arrow115(__unit: None = None, index: Any = index, y: Any = y, xs: Any = xs) -> FSharpList[_T]:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    return reverse(result if is_done else (FSharpList_Cons_305B8EAC(y, result) if ((i + 1) == index) else _arrow115()))


def insert_many_at(index: int, ys: IEnumerable_1[_T], xs: FSharpList[_T]) -> FSharpList[_T]:
    i: int = -1
    is_done: bool = False
    ys_1: FSharpList[_T] = of_seq(ys)

    def folder(acc: FSharpList[_T], x: _T, index: Any = index, ys: Any = ys, xs: Any = xs) -> FSharpList[_T]:
        nonlocal i, is_done
        i = (i + 1) or 0
        if i == index:
            is_done = True
            return FSharpList_Cons_305B8EAC(x, append(ys_1, acc))

        else:
            return FSharpList_Cons_305B8EAC(x, acc)

    result: FSharpList[_T] = fold(folder, FSharpList_get_Empty(), xs)

    def _arrow116(__unit: None = None, index: Any = index, ys: Any = ys, xs: Any = xs) -> FSharpList[_T]:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    return reverse(result if is_done else (append(ys_1, result) if ((i + 1) == index) else _arrow116()))


def remove_at(index: int, xs: FSharpList[_T]) -> FSharpList[_T]:
    i: int = -1
    is_done: bool = False

    def f(_arg: _T | None = None, index: Any = index, xs: Any = xs) -> bool:
        nonlocal i, is_done
        i = (i + 1) or 0
        if i == index:
            is_done = True
            return False

        else:
            return True

    ys: FSharpList[_T] = filter(f, xs)
    if not is_done:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    return ys


def remove_many_at(index: int, count: int, xs: FSharpList[_T]) -> FSharpList[_T]:
    i: int = -1
    status: int = -1

    def f(_arg: _T | None = None, index: Any = index, count: Any = count, xs: Any = xs) -> bool:
        nonlocal i, status
        i = (i + 1) or 0
        if i == index:
            status = 0
            return False

        elif i > index:
            if i < (index + count):
                return False

            else:
                status = 1
                return True

        else:
            return True

    ys: FSharpList[_T] = filter(f, xs)
    status_1: int = (1 if (((i + 1) == (index + count)) if (status == 0) else False) else status) or 0
    if status_1 < 1:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + ("index" if (status_1 < 0) else "count"))

    return ys


def update_at(index: int, y: _T, xs: FSharpList[_T]) -> FSharpList[_T]:
    is_done: bool = False

    def mapping(i: int, x: _T, index: Any = index, y: Any = y, xs: Any = xs) -> _T:
        nonlocal is_done
        if i == index:
            is_done = True
            return y

        else:
            return x

    ys: FSharpList[_T] = map_indexed(mapping, xs)
    if not is_done:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    return ys


__all__ = [
    "FSharpList_reflection",
    "ListEnumerator_1_reflection",
    "FSharpList_get_Empty",
    "FSharpList_Cons_305B8EAC",
    "FSharpList__get_IsEmpty",
    "FSharpList__get_Length",
    "FSharpList__get_Head",
    "FSharpList__get_Tail",
    "FSharpList__get_Item_Z524259A4",
    "index_not_found",
    "empty",
    "cons",
    "singleton",
    "is_empty",
    "length",
    "head",
    "try_head",
    "tail",
    "try_last",
    "last",
    "compare_with",
    "to_array",
    "fold",
    "reverse",
    "fold_back",
    "fold_indexed",
    "fold2",
    "fold_back2",
    "unfold",
    "iterate",
    "iterate2",
    "iterate_indexed",
    "iterate_indexed2",
    "to_seq",
    "of_array_with_tail",
    "of_array",
    "of_seq",
    "concat",
    "scan",
    "scan_back",
    "append",
    "collect",
    "map_indexed",
    "map",
    "indexed",
    "map2",
    "map_indexed2",
    "map3",
    "map_fold",
    "map_fold_back",
    "try_pick",
    "pick",
    "try_find",
    "find",
    "try_find_back",
    "find_back",
    "try_find_index",
    "find_index",
    "try_find_index_back",
    "find_index_back",
    "try_item",
    "item",
    "filter",
    "partition",
    "choose",
    "contains",
    "initialize",
    "replicate",
    "reduce",
    "reduce_back",
    "for_all",
    "for_all2",
    "exists",
    "exists2",
    "unzip",
    "unzip3",
    "zip",
    "zip3",
    "sort_with",
    "sort",
    "sort_by",
    "sort_descending",
    "sort_by_descending",
    "sum",
    "sum_by",
    "max_by",
    "max",
    "min_by",
    "min",
    "average",
    "average_by",
    "permute",
    "chunk_by_size",
    "all_pairs",
    "skip",
    "skip_while",
    "take",
    "take_while",
    "truncate",
    "get_slice",
    "split_at",
    "exactly_one",
    "try_exactly_one",
    "where",
    "pairwise",
    "windowed",
    "split_into",
    "transpose",
    "insert_at",
    "insert_many_at",
    "remove_at",
    "remove_many_at",
    "update_at",
]
