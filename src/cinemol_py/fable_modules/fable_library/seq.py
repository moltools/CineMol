from __future__ import annotations

from collections.abc import Callable
from typing import Any, Generic, TypeVar

from .array_ import chunk_by_size as chunk_by_size_1
from .array_ import fold_back as fold_back_1
from .array_ import fold_back2 as fold_back2_1
from .array_ import map as map_1
from .array_ import map_fold as map_fold_1
from .array_ import map_fold_back as map_fold_back_1
from .array_ import pairwise as pairwise_1
from .array_ import permute as permute_1
from .array_ import reduce_back as reduce_back_1
from .array_ import reverse as reverse_1
from .array_ import scan_back as scan_back_1
from .array_ import singleton as singleton_1
from .array_ import sort_in_place_with
from .array_ import split_into as split_into_1
from .array_ import transpose as transpose_1
from .array_ import try_find_back as try_find_back_1
from .array_ import try_find_index_back as try_find_index_back_1
from .array_ import try_head as try_head_1
from .array_ import try_item as try_item_1
from .array_ import windowed as windowed_1
from .fsharp_core import Operators_NullArg
from .global_ import IGenericAdder_1, IGenericAverager_1, SR_indexOutOfBounds
from .list import FSharpList
from .list import is_empty as is_empty_1
from .list import length as length_1
from .list import of_array as of_array_1
from .list import of_seq as of_seq_1
from .list import to_array as to_array_1
from .list import try_head as try_head_2
from .list import try_item as try_item_2
from .option import some
from .option import value as value_1
from .reflection import TypeInfo, class_type
from .types import Array, to_string
from .util import (
    IComparer_1,
    IDisposable,
    IEnumerable_1,
    IEnumerator,
    IEqualityComparer_1,
    clear,
    equals,
    get_enumerator,
    ignore,
    is_array_like,
    is_disposable,
    lock,
    to_enumerable,
    to_iterator,
)
from .util import dispose as dispose_2


_T = TypeVar("_T")

_U = TypeVar("_U")

_STATE = TypeVar("_STATE")

__B = TypeVar("__B")

__A = TypeVar("__A")

_T1 = TypeVar("_T1")

_T2 = TypeVar("_T2")

_T3 = TypeVar("_T3")

_RESULT = TypeVar("_RESULT")

_COLLECTION = TypeVar("_COLLECTION")

SR_enumerationAlreadyFinished: str = "Enumeration already finished."

SR_enumerationNotStarted: str = "Enumeration has not started. Call MoveNext."

SR_inputSequenceEmpty: str = "The input sequence was empty."

SR_inputSequenceTooLong: str = "The input sequence contains more than one element."

SR_keyNotFoundAlt: str = "An index satisfying the predicate was not found in the collection."

SR_notEnoughElements: str = "The input sequence has an insufficient number of elements."

SR_resetNotSupported: str = "Reset is not supported on this enumerator."


def Enumerator_noReset(__unit: None = None) -> Any:
    raise Exception(SR_resetNotSupported)


def Enumerator_notStarted(__unit: None = None) -> Any:
    raise Exception(SR_enumerationNotStarted)


def Enumerator_alreadyFinished(__unit: None = None) -> Any:
    raise Exception(SR_enumerationAlreadyFinished)


def _expr136(gen0: TypeInfo) -> TypeInfo:
    return class_type("SeqModule.Enumerator.Seq", [gen0], Enumerator_Seq)


class Enumerator_Seq(Generic[_T]):
    def __init__(self, f: Callable[[], IEnumerator[_T]]) -> None:
        self.f: Callable[[], IEnumerator[_T]] = f

    def __str__(self, __unit: None = None) -> str:
        xs: Enumerator_Seq[_T] = self
        i: int = 0
        str_1: str = "seq ["
        with get_enumerator(xs) as e:
            while e.System_Collections_IEnumerator_MoveNext() if (i < 4) else False:
                if i > 0:
                    str_1 = str_1 + "; "

                str_1 = str_1 + to_string(e.System_Collections_Generic_IEnumerator_1_get_Current())
                i = (i + 1) or 0
            if i == 4:
                str_1 = str_1 + "; ..."

            return str_1 + "]"

    def GetEnumerator(self, __unit: None = None) -> IEnumerator[_T]:
        x: Enumerator_Seq[_T] = self
        return x.f(None)

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None = None) -> IEnumerator[Any]:
        x: Enumerator_Seq[_T] = self
        return x.f(None)


Enumerator_Seq_reflection = _expr136


def Enumerator_Seq__ctor_673A07F2(f: Callable[[], IEnumerator[_T]]) -> Enumerator_Seq[_T]:
    return Enumerator_Seq(f)


def _expr137(gen0: TypeInfo) -> TypeInfo:
    return class_type("SeqModule.Enumerator.FromFunctions`1", [gen0], Enumerator_FromFunctions_1)


class Enumerator_FromFunctions_1(IDisposable, Generic[_T]):
    def __init__(self, current: Callable[[], _T], next_1: Callable[[], bool], dispose: Callable[[], None]) -> None:
        self.current: Callable[[], _T] = current
        self.next: Callable[[], bool] = next_1
        self.dispose: Callable[[], None] = dispose

    def System_Collections_Generic_IEnumerator_1_get_Current(self, __unit: None = None) -> _T:
        _: Enumerator_FromFunctions_1[_T] = self
        return _.current(None)

    def System_Collections_IEnumerator_get_Current(self, __unit: None = None) -> Any:
        _: Enumerator_FromFunctions_1[_T] = self
        return _.current(None)

    def System_Collections_IEnumerator_MoveNext(self, __unit: None = None) -> bool:
        _: Enumerator_FromFunctions_1[_T] = self
        return _.next(None)

    def System_Collections_IEnumerator_Reset(self, __unit: None = None) -> None:
        Enumerator_noReset()

    def Dispose(self, __unit: None = None) -> None:
        _: Enumerator_FromFunctions_1[_T] = self
        _.dispose(None)


Enumerator_FromFunctions_1_reflection = _expr137


def Enumerator_FromFunctions_1__ctor_58C54629(
    current: Callable[[], _T], next_1: Callable[[], bool], dispose: Callable[[], None]
) -> Enumerator_FromFunctions_1[_T]:
    return Enumerator_FromFunctions_1(current, next_1, dispose)


def Enumerator_cast(e: IEnumerator[_T]) -> IEnumerator[_T]:
    def current(__unit: None = None, e: Any = e) -> _T:
        return e.System_Collections_Generic_IEnumerator_1_get_Current()

    def next_1(__unit: None = None, e: Any = e) -> bool:
        return e.System_Collections_IEnumerator_MoveNext()

    def dispose(__unit: None = None, e: Any = e) -> None:
        dispose_2(e)

    return Enumerator_FromFunctions_1__ctor_58C54629(current, next_1, dispose)


def Enumerator_concat(sources: IEnumerable_1[Any]) -> IEnumerator[Any]:
    outer_opt: IEnumerator[_U] | None = None
    inner_opt: IEnumerator[_T] | None = None
    started: bool = False
    finished: bool = False
    curr: _T | None = None

    def finish(__unit: None = None, sources: Any = sources) -> None:
        nonlocal finished, inner_opt, outer_opt
        finished = True
        if inner_opt is not None:
            inner: IEnumerator[_T] = inner_opt
            try:
                dispose_2(inner)

            finally:
                inner_opt = None

        if outer_opt is not None:
            outer: IEnumerator[_U] = outer_opt
            try:
                dispose_2(outer)

            finally:
                outer_opt = None

    def current(__unit: None = None, sources: Any = sources) -> _T:
        if not started:
            Enumerator_notStarted()

        elif finished:
            Enumerator_alreadyFinished()

        if curr is not None:
            return value_1(curr)

        else:
            return Enumerator_alreadyFinished()

    def next_1(__unit: None = None, sources: Any = sources) -> bool:
        nonlocal started
        if not started:
            started = True

        if finished:
            return False

        else:
            res: bool | None = None
            while res is None:
                nonlocal curr, inner_opt, outer_opt
                outer_opt_1: IEnumerator[_U] | None = outer_opt
                inner_opt_1: IEnumerator[_T] | None = inner_opt
                if outer_opt_1 is not None:
                    if inner_opt_1 is not None:
                        inner_1: IEnumerator[_T] = inner_opt_1
                        if inner_1.System_Collections_IEnumerator_MoveNext():
                            curr = some(inner_1.System_Collections_Generic_IEnumerator_1_get_Current())
                            res = True

                        else:
                            try:
                                dispose_2(inner_1)

                            finally:
                                inner_opt = None

                    else:
                        outer_1: IEnumerator[_U] = outer_opt_1
                        if outer_1.System_Collections_IEnumerator_MoveNext():
                            ie: _U = outer_1.System_Collections_Generic_IEnumerator_1_get_Current()

                            def _arrow138(__unit: None = None) -> IEnumerator[_T]:
                                copy_of_struct: _U = ie
                                return get_enumerator(copy_of_struct)

                            inner_opt = _arrow138()

                        else:
                            finish(None)
                            res = False

                else:
                    outer_opt = get_enumerator(sources)

            return value_1(res)

    def dispose(__unit: None = None, sources: Any = sources) -> None:
        if not finished:
            finish(None)

    return Enumerator_FromFunctions_1__ctor_58C54629(current, next_1, dispose)


def Enumerator_enumerateThenFinally(f: Callable[[], None], e: IEnumerator[_T]) -> IEnumerator[_T]:
    def current(__unit: None = None, f: Any = f, e: Any = e) -> _T:
        return e.System_Collections_Generic_IEnumerator_1_get_Current()

    def next_1(__unit: None = None, f: Any = f, e: Any = e) -> bool:
        return e.System_Collections_IEnumerator_MoveNext()

    def dispose(__unit: None = None, f: Any = f, e: Any = e) -> None:
        try:
            dispose_2(e)

        finally:
            f(None)

    return Enumerator_FromFunctions_1__ctor_58C54629(current, next_1, dispose)


def Enumerator_generateWhileSome(
    openf: Callable[[], _T], compute: Callable[[_T], _U | None], closef: Callable[[_T], None]
) -> IEnumerator[_U]:
    started: bool = False
    curr: _U | None = None
    state: _T | None = some(openf(None))

    def dispose(__unit: None = None, openf: Any = openf, compute: Any = compute, closef: Any = closef) -> None:
        nonlocal state
        if state is not None:
            x_1: _T = value_1(state)
            try:
                closef(x_1)

            finally:
                state = None

    def finish(__unit: None = None, openf: Any = openf, compute: Any = compute, closef: Any = closef) -> None:
        nonlocal curr
        try:
            dispose(None)

        finally:
            curr = None

    def current(__unit: None = None, openf: Any = openf, compute: Any = compute, closef: Any = closef) -> _U:
        if not started:
            Enumerator_notStarted()

        if curr is not None:
            return value_1(curr)

        else:
            return Enumerator_alreadyFinished()

    def next_1(__unit: None = None, openf: Any = openf, compute: Any = compute, closef: Any = closef) -> bool:
        nonlocal started, curr
        if not started:
            started = True

        if state is not None:
            s: _T = value_1(state)
            match_value_1: _U | None
            try:
                match_value_1 = compute(s)

            except Exception as match_value:
                finish(None)
                raise match_value

            if match_value_1 is not None:
                curr = match_value_1
                return True

            else:
                finish(None)
                return False

        else:
            return False

    return Enumerator_FromFunctions_1__ctor_58C54629(current, next_1, dispose)


def Enumerator_unfold(f: Callable[[_STATE], tuple[_T, _STATE] | None], state: _STATE) -> IEnumerator[_T]:
    curr: tuple[_T, _STATE] | None = None
    acc: _STATE = state

    def current(__unit: None = None, f: Any = f, state: Any = state) -> _T:
        if curr is not None:
            x: _T = curr[0]
            st: _STATE = curr[1]
            return x

        else:
            return Enumerator_notStarted()

    def next_1(__unit: None = None, f: Any = f, state: Any = state) -> bool:
        nonlocal curr, acc
        curr = f(acc)
        if curr is not None:
            x_1: _T = curr[0]
            st_1: _STATE = curr[1]
            acc = st_1
            return True

        else:
            return False

    def dispose(__unit: None = None, f: Any = f, state: Any = state) -> None:
        pass

    return Enumerator_FromFunctions_1__ctor_58C54629(current, next_1, dispose)


def index_not_found(__unit: None = None) -> Any:
    raise Exception(SR_keyNotFoundAlt)


def check_non_null(arg_name: str, arg: Any) -> None:
    if arg is None:
        Operators_NullArg(arg_name)


def mk_seq(f: Callable[[], IEnumerator[_T]]) -> IEnumerable_1[_T]:
    return Enumerator_Seq__ctor_673A07F2(f)


def of_seq(xs: IEnumerable_1[_T]) -> IEnumerator[_T]:
    check_non_null("source", xs)
    return get_enumerator(xs)


def delay(generator: Callable[[], IEnumerable_1[_T]]) -> IEnumerable_1[_T]:
    def _arrow139(__unit: None = None, generator: Any = generator) -> IEnumerator[_T]:
        return get_enumerator(generator(None))

    return mk_seq(_arrow139)


def concat(sources: IEnumerable_1[Any]) -> IEnumerable_1[Any]:
    def _arrow140(__unit: None = None, sources: Any = sources) -> IEnumerator[_T]:
        return Enumerator_concat(sources)

    return mk_seq(_arrow140)


def unfold(generator: Callable[[_STATE], tuple[_T, _STATE] | None], state: _STATE) -> IEnumerable_1[_T]:
    def _arrow141(__unit: None = None, generator: Any = generator, state: Any = state) -> IEnumerator[_T]:
        return Enumerator_unfold(generator, state)

    return mk_seq(_arrow141)


def empty(__unit: None = None) -> IEnumerable_1[Any]:
    def _arrow142(__unit: None = None) -> IEnumerable_1[_T]:
        return [0] * 0

    return delay(_arrow142)


def singleton(x: _T | None = None) -> IEnumerable_1[_T]:
    def _arrow143(__unit: None = None, x: Any = x) -> IEnumerable_1[_T]:
        return singleton_1(x, None)

    return delay(_arrow143)


def of_array(arr: Array[_T]) -> IEnumerable_1[_T]:
    return arr


def to_array(xs: IEnumerable_1[_T]) -> Array[_T]:
    if isinstance(xs, FSharpList):
        return to_array_1(xs)

    else:
        return list(xs)


def of_list(xs: FSharpList[_T]) -> IEnumerable_1[_T]:
    return xs


def to_list(xs: IEnumerable_1[_T]) -> FSharpList[_T]:
    if is_array_like(xs):
        return of_array_1(xs)

    elif isinstance(xs, FSharpList):
        return xs

    else:
        return of_seq_1(xs)


def generate(
    create: Callable[[], __A], compute: Callable[[__A], __B | None], dispose: Callable[[__A], None]
) -> IEnumerable_1[__B]:
    def _arrow144(
        __unit: None = None, create: Any = create, compute: Any = compute, dispose: Any = dispose
    ) -> IEnumerator[__B]:
        return Enumerator_generateWhileSome(create, compute, dispose)

    return mk_seq(_arrow144)


def generate_indexed(
    create: Callable[[], __A], compute: Callable[[int, __A], __B | None], dispose: Callable[[__A], None]
) -> IEnumerable_1[__B]:
    def _arrow146(
        __unit: None = None, create: Any = create, compute: Any = compute, dispose: Any = dispose
    ) -> IEnumerator[__B]:
        i: int = -1

        def _arrow145(x: __A | None = None) -> __B | None:
            nonlocal i
            i = (i + 1) or 0
            return compute(i, x)

        return Enumerator_generateWhileSome(create, _arrow145, dispose)

    return mk_seq(_arrow146)


def append(xs: IEnumerable_1[_T], ys: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    return concat(to_enumerable([xs, ys]))


def cast(xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow147(__unit: None = None, xs: Any = xs) -> IEnumerator[_T]:
        check_non_null("source", xs)
        return Enumerator_cast(get_enumerator(xs))

    return mk_seq(_arrow147)


def choose(chooser: Callable[[_T], _U | None], xs: IEnumerable_1[_T]) -> IEnumerable_1[_U]:
    def _arrow148(__unit: None = None, chooser: Any = chooser, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow149(e: IEnumerator[_T], chooser: Any = chooser, xs: Any = xs) -> _U | None:
        curr: _U | None = None
        while e.System_Collections_IEnumerator_MoveNext() if (curr is None) else False:
            curr = chooser(e.System_Collections_Generic_IEnumerator_1_get_Current())
        return curr

    def _arrow150(e_1: IEnumerator[_T], chooser: Any = chooser, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate(_arrow148, _arrow149, _arrow150)


def compare_with(comparer: Callable[[_T, _T], int], xs: IEnumerable_1[_T], ys: IEnumerable_1[_T]) -> int:
    with of_seq(xs) as e1:
        with of_seq(ys) as e2:
            c: int = 0
            b1: bool = e1.System_Collections_IEnumerator_MoveNext()
            b2: bool = e2.System_Collections_IEnumerator_MoveNext()
            while b2 if (b1 if (c == 0) else False) else False:
                c = (
                    comparer(
                        e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                        e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                    )
                    or 0
                )
                if c == 0:
                    b1 = e1.System_Collections_IEnumerator_MoveNext()
                    b2 = e2.System_Collections_IEnumerator_MoveNext()

            if c != 0:
                return c

            elif b1:
                return 1

            elif b2:
                return -1

            else:
                return 0


def contains(value: _T, xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> bool:
    with of_seq(xs) as e:
        found: bool = False
        while e.System_Collections_IEnumerator_MoveNext() if (not found) else False:
            found = comparer.Equals(value, e.System_Collections_Generic_IEnumerator_1_get_Current())
        return found


def enumerate_from_functions(
    create: Callable[[], __A], move_next: Callable[[__A], bool], current: Callable[[__A], __B]
) -> IEnumerable_1[__B]:
    def _arrow151(
        x: __A | None = None, create: Any = create, move_next: Any = move_next, current: Any = current
    ) -> __B | None:
        return some(current(x)) if move_next(x) else None

    def _arrow152(
        x_1: __A | None = None, create: Any = create, move_next: Any = move_next, current: Any = current
    ) -> None:
        match_value: Any = x_1
        if is_disposable(match_value):
            dispose_2(match_value)

    return generate(create, _arrow151, _arrow152)


def enumerate_then_finally(source: IEnumerable_1[_T], compensation: Callable[[], None]) -> IEnumerable_1[_T]:
    compensation_1: Callable[[], None] = compensation

    def _arrow153(__unit: None = None, source: Any = source, compensation: Any = compensation) -> IEnumerator[_T]:
        try:
            return Enumerator_enumerateThenFinally(compensation_1, of_seq(source))

        except Exception as match_value:
            compensation_1(None)
            raise match_value

    return mk_seq(_arrow153)


def enumerate_using(resource: _T, source: Callable[[_T], __A]) -> IEnumerable_1[Any]:
    def compensation(__unit: None = None, resource: Any = resource, source: Any = source) -> None:
        if equals(resource, None):
            pass

        else:
            copy_of_struct: _T = resource
            dispose_2(copy_of_struct)

    def _arrow154(__unit: None = None, resource: Any = resource, source: Any = source) -> IEnumerator[_U]:
        try:
            return Enumerator_enumerateThenFinally(compensation, of_seq(source(resource)))

        except Exception as match_value_1:
            compensation(None)
            raise match_value_1

    return mk_seq(_arrow154)


def enumerate_while(guard: Callable[[], bool], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow155(i: int, guard: Any = guard, xs: Any = xs) -> tuple[IEnumerable_1[_T], int] | None:
        return ((xs, i + 1)) if guard(None) else None

    return concat(unfold(_arrow155, 0))


def filter(f: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def chooser(x: _T | None = None, f: Any = f, xs: Any = xs) -> _T | None:
        if f(x):
            return some(x)

        else:
            return None

    return choose(chooser, xs)


def exists(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> bool:
    with of_seq(xs) as e:
        found: bool = False
        while e.System_Collections_IEnumerator_MoveNext() if (not found) else False:
            found = predicate(e.System_Collections_Generic_IEnumerator_1_get_Current())
        return found


def exists2(predicate: Callable[[_T1, _T2], bool], xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]) -> bool:
    with of_seq(xs) as e1:
        with of_seq(ys) as e2:
            found: bool = False
            while (
                e2.System_Collections_IEnumerator_MoveNext()
                if (e1.System_Collections_IEnumerator_MoveNext() if (not found) else False)
                else False
            ):
                found = predicate(
                    e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                )
            return found


def exactly_one(xs: IEnumerable_1[_T]) -> _T:
    with of_seq(xs) as e:
        if e.System_Collections_IEnumerator_MoveNext():
            v: _T = e.System_Collections_Generic_IEnumerator_1_get_Current()
            if e.System_Collections_IEnumerator_MoveNext():
                raise Exception((SR_inputSequenceTooLong + "\\nParameter name: ") + "source")

            else:
                return v

        else:
            raise Exception((SR_inputSequenceEmpty + "\\nParameter name: ") + "source")


def try_exactly_one(xs: IEnumerable_1[_T]) -> _T | None:
    with of_seq(xs) as e:
        if e.System_Collections_IEnumerator_MoveNext():
            v: _T = e.System_Collections_Generic_IEnumerator_1_get_Current()
            if e.System_Collections_IEnumerator_MoveNext():
                return None

            else:
                return some(v)

        else:
            return None


def try_find(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> _T | None:
    with of_seq(xs) as e:
        res: _T | None = None
        while e.System_Collections_IEnumerator_MoveNext() if (res is None) else False:
            c: _T = e.System_Collections_Generic_IEnumerator_1_get_Current()
            if predicate(c):
                res = some(c)

        return res


def find(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> _T:
    match_value: _T | None = try_find(predicate, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def try_find_back(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> _T | None:
    return try_find_back_1(predicate, to_array(xs))


def find_back(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> _T:
    match_value: _T | None = try_find_back(predicate, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def try_find_index(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> int | None:
    with of_seq(xs) as e:

        def loop(i_mut: int, predicate: Any = predicate, xs: Any = xs) -> int | None:
            while True:
                (i,) = (i_mut,)
                if e.System_Collections_IEnumerator_MoveNext():
                    if predicate(e.System_Collections_Generic_IEnumerator_1_get_Current()):
                        return i

                    else:
                        i_mut = i + 1
                        continue

                else:
                    return None

                break

        return loop(0)


def find_index(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> int:
    match_value: int | None = try_find_index(predicate, xs)
    if match_value is None:
        index_not_found()
        return -1

    else:
        return match_value


def try_find_index_back(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> int | None:
    return try_find_index_back_1(predicate, to_array(xs))


def find_index_back(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> int:
    match_value: int | None = try_find_index_back(predicate, xs)
    if match_value is None:
        index_not_found()
        return -1

    else:
        return match_value


def fold(folder: Callable[[_STATE, _T], _STATE], state: _STATE, xs: IEnumerable_1[_T]) -> _STATE:
    with of_seq(xs) as e:
        acc: _STATE = state
        while e.System_Collections_IEnumerator_MoveNext():
            acc = folder(acc, e.System_Collections_Generic_IEnumerator_1_get_Current())
        return acc


def fold_back(folder: Callable[[_T, Any], Any], xs: IEnumerable_1[_T], state: Any = None) -> Any:
    return fold_back_1(folder, to_array(xs), state)


def fold2(
    folder: Callable[[_STATE, _T1, _T2], _STATE], state: _STATE, xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]
) -> _STATE:
    with of_seq(xs) as e1:
        with of_seq(ys) as e2:
            acc: _STATE = state
            while (
                e2.System_Collections_IEnumerator_MoveNext() if e1.System_Collections_IEnumerator_MoveNext() else False
            ):
                acc = folder(
                    acc,
                    e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                )
            return acc


def fold_back2(
    folder: Callable[[_T1, _T2, _STATE], _STATE], xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2], state: _STATE
) -> _STATE:
    return fold_back2_1(folder, to_array(xs), to_array(ys), state)


def for_all(predicate: Callable[[__A], bool], xs: IEnumerable_1[__A]) -> bool:
    def _arrow156(x: __A | None = None, predicate: Any = predicate, xs: Any = xs) -> bool:
        return not predicate(x)

    return not exists(_arrow156, xs)


def for_all2(predicate: Callable[[__A, __B], bool], xs: IEnumerable_1[__A], ys: IEnumerable_1[__B]) -> bool:
    def _arrow157(x: __A, y: __B, predicate: Any = predicate, xs: Any = xs, ys: Any = ys) -> bool:
        return not predicate(x, y)

    return not exists2(_arrow157, xs, ys)


def try_head(xs: IEnumerable_1[_T]) -> _T | None:
    if is_array_like(xs):
        return try_head_1(xs)

    elif isinstance(xs, FSharpList):
        return try_head_2(xs)

    else:
        with of_seq(xs) as e:
            if e.System_Collections_IEnumerator_MoveNext():
                return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

            else:
                return None


def head(xs: IEnumerable_1[_T]) -> _T:
    match_value: _T | None = try_head(xs)
    if match_value is None:
        raise Exception((SR_inputSequenceEmpty + "\\nParameter name: ") + "source")

    else:
        return value_1(match_value)


def initialize(count: int, f: Callable[[int], __A]) -> IEnumerable_1[__A]:
    def _arrow158(i: int, count: Any = count, f: Any = f) -> tuple[__A, int] | None:
        return ((f(i), i + 1)) if (i < count) else None

    return unfold(_arrow158, 0)


def initialize_infinite(f: Callable[[int], __A]) -> IEnumerable_1[__A]:
    return initialize(2147483647, f)


def is_empty(xs: IEnumerable_1[Any]) -> bool:
    if is_array_like(xs):
        return len(xs) == 0

    elif isinstance(xs, FSharpList):
        return is_empty_1(xs)

    else:
        with of_seq(xs) as e:
            return not e.System_Collections_IEnumerator_MoveNext()


def try_item(index: int, xs: IEnumerable_1[_T]) -> _T | None:
    if is_array_like(xs):
        return try_item_1(index, xs)

    elif isinstance(xs, FSharpList):
        return try_item_2(index, xs)

    else:
        with of_seq(xs) as e:

            def loop(index_1_mut: int, index: Any = index, xs: Any = xs) -> _T | None:
                while True:
                    (index_1,) = (index_1_mut,)
                    if not e.System_Collections_IEnumerator_MoveNext():
                        return None

                    elif index_1 == 0:
                        return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

                    else:
                        index_1_mut = index_1 - 1
                        continue

                    break

            return loop(index)


def item(index: int, xs: IEnumerable_1[_T]) -> _T:
    match_value: _T | None = try_item(index, xs)
    if match_value is None:
        raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "index")

    else:
        return value_1(match_value)


def iterate(action: Callable[[__A], None], xs: IEnumerable_1[__A]) -> None:
    def _arrow159(unit_var: None, x: __A, action: Any = action, xs: Any = xs) -> None:
        action(x)

    fold(_arrow159, None, xs)


def iterate2(action: Callable[[__A, __B], None], xs: IEnumerable_1[__A], ys: IEnumerable_1[__B]) -> None:
    def _arrow160(unit_var: None, x: __A, y: __B, action: Any = action, xs: Any = xs, ys: Any = ys) -> None:
        action(x, y)

    fold2(_arrow160, None, xs, ys)


def iterate_indexed(action: Callable[[int, __A], None], xs: IEnumerable_1[__A]) -> None:
    def _arrow161(i: int, x: __A, action: Any = action, xs: Any = xs) -> int:
        action(i, x)
        return i + 1

    ignore(fold(_arrow161, 0, xs))


def iterate_indexed2(action: Callable[[int, __A, __B], None], xs: IEnumerable_1[__A], ys: IEnumerable_1[__B]) -> None:
    def _arrow162(i: int, x: __A, y: __B, action: Any = action, xs: Any = xs, ys: Any = ys) -> int:
        action(i, x, y)
        return i + 1

    ignore(fold2(_arrow162, 0, xs, ys))


def try_last(xs: IEnumerable_1[_T]) -> _T | None:
    with of_seq(xs) as e:

        def loop(acc_mut: _T | None = None, xs: Any = xs) -> _T | None:
            while True:
                (acc,) = (acc_mut,)
                if not e.System_Collections_IEnumerator_MoveNext():
                    return acc

                else:
                    acc_mut = e.System_Collections_Generic_IEnumerator_1_get_Current()
                    continue

                break

        if e.System_Collections_IEnumerator_MoveNext():
            return some(loop(e.System_Collections_Generic_IEnumerator_1_get_Current()))

        else:
            return None


def last(xs: IEnumerable_1[_T]) -> _T:
    match_value: _T | None = try_last(xs)
    if match_value is None:
        raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "source")

    else:
        return value_1(match_value)


def length(xs: IEnumerable_1[Any]) -> int:
    if is_array_like(xs):
        return len(xs)

    elif isinstance(xs, FSharpList):
        return length_1(xs)

    else:
        with of_seq(xs) as e:
            count: int = 0
            while e.System_Collections_IEnumerator_MoveNext():
                count = (count + 1) or 0
            return count


def map(mapping: Callable[[_T], _U], xs: IEnumerable_1[_T]) -> IEnumerable_1[_U]:
    def _arrow163(__unit: None = None, mapping: Any = mapping, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow164(e: IEnumerator[_T], mapping: Any = mapping, xs: Any = xs) -> _U | None:
        return (
            some(mapping(e.System_Collections_Generic_IEnumerator_1_get_Current()))
            if e.System_Collections_IEnumerator_MoveNext()
            else None
        )

    def _arrow165(e_1: IEnumerator[_T], mapping: Any = mapping, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate(_arrow163, _arrow164, _arrow165)


def map_indexed(mapping: Callable[[int, _T], _U], xs: IEnumerable_1[_T]) -> IEnumerable_1[_U]:
    def _arrow166(__unit: None = None, mapping: Any = mapping, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow167(i: int, e: IEnumerator[_T], mapping: Any = mapping, xs: Any = xs) -> _U | None:
        return (
            some(mapping(i, e.System_Collections_Generic_IEnumerator_1_get_Current()))
            if e.System_Collections_IEnumerator_MoveNext()
            else None
        )

    def _arrow168(e_1: IEnumerator[_T], mapping: Any = mapping, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow166, _arrow167, _arrow168)


def indexed(xs: IEnumerable_1[_T]) -> IEnumerable_1[tuple[int, _T]]:
    def mapping(i: int, x: _T, xs: Any = xs) -> tuple[int, _T]:
        return (i, x)

    return map_indexed(mapping, xs)


def map2(mapping: Callable[[_T1, _T2], _U], xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]) -> IEnumerable_1[_U]:
    def _arrow169(
        __unit: None = None, mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> tuple[IEnumerator[_T1], IEnumerator[_T2]]:
        return (of_seq(xs), of_seq(ys))

    def _arrow170(
        tupled_arg: tuple[IEnumerator[_T1], IEnumerator[_T2]], mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> _U | None:
        e1: IEnumerator[_T1] = tupled_arg[0]
        e2: IEnumerator[_T2] = tupled_arg[1]
        return (
            some(
                mapping(
                    e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                )
            )
            if (e2.System_Collections_IEnumerator_MoveNext() if e1.System_Collections_IEnumerator_MoveNext() else False)
            else None
        )

    def _arrow171(
        tupled_arg_1: tuple[IEnumerator[_T1], IEnumerator[_T2]], mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> None:
        try:
            dispose_2(tupled_arg_1[0])

        finally:
            dispose_2(tupled_arg_1[1])

    return generate(_arrow169, _arrow170, _arrow171)


def map_indexed2(
    mapping: Callable[[int, _T1, _T2], _U], xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]
) -> IEnumerable_1[_U]:
    def _arrow172(
        __unit: None = None, mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> tuple[IEnumerator[_T1], IEnumerator[_T2]]:
        return (of_seq(xs), of_seq(ys))

    def _arrow173(
        i: int,
        tupled_arg: tuple[IEnumerator[_T1], IEnumerator[_T2]],
        mapping: Any = mapping,
        xs: Any = xs,
        ys: Any = ys,
    ) -> _U | None:
        e1: IEnumerator[_T1] = tupled_arg[0]
        e2: IEnumerator[_T2] = tupled_arg[1]
        return (
            some(
                mapping(
                    i,
                    e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                )
            )
            if (e2.System_Collections_IEnumerator_MoveNext() if e1.System_Collections_IEnumerator_MoveNext() else False)
            else None
        )

    def _arrow174(
        tupled_arg_1: tuple[IEnumerator[_T1], IEnumerator[_T2]], mapping: Any = mapping, xs: Any = xs, ys: Any = ys
    ) -> None:
        try:
            dispose_2(tupled_arg_1[0])

        finally:
            dispose_2(tupled_arg_1[1])

    return generate_indexed(_arrow172, _arrow173, _arrow174)


def map3(
    mapping: Callable[[_T1, _T2, _T3], _U], xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2], zs: IEnumerable_1[_T3]
) -> IEnumerable_1[_U]:
    def _arrow175(
        __unit: None = None, mapping: Any = mapping, xs: Any = xs, ys: Any = ys, zs: Any = zs
    ) -> tuple[IEnumerator[_T1], IEnumerator[_T2], IEnumerator[_T3]]:
        return (of_seq(xs), of_seq(ys), of_seq(zs))

    def _arrow176(
        tupled_arg: tuple[IEnumerator[_T1], IEnumerator[_T2], IEnumerator[_T3]],
        mapping: Any = mapping,
        xs: Any = xs,
        ys: Any = ys,
        zs: Any = zs,
    ) -> _U | None:
        e1: IEnumerator[_T1] = tupled_arg[0]
        e2: IEnumerator[_T2] = tupled_arg[1]
        e3: IEnumerator[_T3] = tupled_arg[2]
        return (
            some(
                mapping(
                    e1.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e2.System_Collections_Generic_IEnumerator_1_get_Current(),
                    e3.System_Collections_Generic_IEnumerator_1_get_Current(),
                )
            )
            if (
                e3.System_Collections_IEnumerator_MoveNext()
                if (
                    e2.System_Collections_IEnumerator_MoveNext()
                    if e1.System_Collections_IEnumerator_MoveNext()
                    else False
                )
                else False
            )
            else None
        )

    def _arrow177(
        tupled_arg_1: tuple[IEnumerator[_T1], IEnumerator[_T2], IEnumerator[_T3]],
        mapping: Any = mapping,
        xs: Any = xs,
        ys: Any = ys,
        zs: Any = zs,
    ) -> None:
        try:
            dispose_2(tupled_arg_1[0])

        finally:
            try:
                dispose_2(tupled_arg_1[1])

            finally:
                dispose_2(tupled_arg_1[2])

    return generate(_arrow175, _arrow176, _arrow177)


def read_only(xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    check_non_null("source", xs)

    def _arrow178(x: _T | None = None, xs: Any = xs) -> _T | None:
        return x

    return map(_arrow178, xs)


def _expr179(gen0: TypeInfo) -> TypeInfo:
    return class_type("SeqModule.CachedSeq`1", [gen0], CachedSeq_1)


class CachedSeq_1(IDisposable, Generic[_T]):
    def __init__(self, cleanup: Callable[[], None], res: IEnumerable_1[Any]) -> None:
        self.cleanup: Callable[[], None] = cleanup
        self.res: IEnumerable_1[_T] = res

    def Dispose(self, __unit: None = None) -> None:
        _: CachedSeq_1[_T] = self
        _.cleanup(None)

    def GetEnumerator(self, __unit: None = None) -> IEnumerator[_T]:
        _: CachedSeq_1[_T] = self
        return get_enumerator(_.res)

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None = None) -> IEnumerator[Any]:
        _: CachedSeq_1[_T] = self
        return get_enumerator(_.res)


CachedSeq_1_reflection = _expr179


def CachedSeq_1__ctor_Z7A8347D4(cleanup: Callable[[], None], res: IEnumerable_1[Any]) -> CachedSeq_1[_T]:
    return CachedSeq_1(cleanup, res)


def CachedSeq_1__Clear(_: CachedSeq_1[Any]) -> None:
    _.cleanup(None)


def cache(source: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    check_non_null("source", source)
    prefix: Array[_T] = []
    enumerator_r: (IEnumerator[_T] | None) | None = None

    def cleanup(__unit: None = None, source: Any = source) -> None:
        def action_1(__unit: None = None) -> None:
            nonlocal enumerator_r
            clear(prefix)
            (pattern_matching_result, e) = (None, None)
            if enumerator_r is not None:
                if value_1(enumerator_r) is not None:
                    pattern_matching_result = 0
                    e = value_1(enumerator_r)

                else:
                    pattern_matching_result = 1

            else:
                pattern_matching_result = 1

            if pattern_matching_result == 0:
                dispose_2(e)

            enumerator_r = None

        lock(prefix, action_1)

    def _arrow180(i_1: int, source: Any = source) -> tuple[_T, int] | None:
        def action(__unit: None = None) -> tuple[_T, int] | None:
            nonlocal enumerator_r
            if i_1 < len(prefix):
                return (prefix[i_1], i_1 + 1)

            else:
                if i_1 >= len(prefix):
                    opt_enumerator_2: IEnumerator[_T] | None
                    if enumerator_r is not None:
                        opt_enumerator_2 = value_1(enumerator_r)

                    else:
                        opt_enumerator: IEnumerator[_T] | None = get_enumerator(source)
                        enumerator_r = some(opt_enumerator)
                        opt_enumerator_2 = opt_enumerator

                    if opt_enumerator_2 is None:
                        pass

                    else:
                        enumerator: IEnumerator[_T] = opt_enumerator_2
                        if enumerator.System_Collections_IEnumerator_MoveNext():
                            (prefix.append(enumerator.System_Collections_Generic_IEnumerator_1_get_Current()))

                        else:
                            dispose_2(enumerator)
                            enumerator_r = some(None)

                if i_1 < len(prefix):
                    return (prefix[i_1], i_1 + 1)

                else:
                    return None

        return lock(prefix, action)

    return CachedSeq_1__ctor_Z7A8347D4(cleanup, unfold(_arrow180, 0))


def all_pairs(xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]) -> IEnumerable_1[tuple[_T1, _T2]]:
    ys_cache: IEnumerable_1[_T2] = cache(ys)

    def _arrow181(__unit: None = None, xs: Any = xs, ys: Any = ys) -> IEnumerable_1[tuple[_T1, _T2]]:
        def mapping_1(x: _T1 | None = None) -> IEnumerable_1[tuple[_T1, _T2]]:
            def mapping(y: _T2 | None = None, x: Any = x) -> tuple[_T1, _T2]:
                return (x, y)

            return map(mapping, ys_cache)

        return concat(map(mapping_1, xs))

    return delay(_arrow181)


def map_fold(
    mapping: Callable[[_STATE, _T], tuple[_RESULT, _STATE]], state: _STATE, xs: IEnumerable_1[_T]
) -> tuple[IEnumerable_1[_RESULT], _STATE]:
    pattern_input: tuple[Array[_RESULT], _STATE] = map_fold_1(mapping, state, to_array(xs), None)
    return (read_only(pattern_input[0]), pattern_input[1])


def map_fold_back(
    mapping: Callable[[_T, _STATE], tuple[_RESULT, _STATE]], xs: IEnumerable_1[_T], state: _STATE
) -> tuple[IEnumerable_1[_RESULT], _STATE]:
    pattern_input: tuple[Array[_RESULT], _STATE] = map_fold_back_1(mapping, to_array(xs), state, None)
    return (read_only(pattern_input[0]), pattern_input[1])


def try_pick(chooser: Callable[[_T], __A | None], xs: IEnumerable_1[_T]) -> __A | None:
    with of_seq(xs) as e:
        res: __A | None = None
        while e.System_Collections_IEnumerator_MoveNext() if (res is None) else False:
            res = chooser(e.System_Collections_Generic_IEnumerator_1_get_Current())
        return res


def pick(chooser: Callable[[_T], __A | None], xs: IEnumerable_1[_T]) -> __A:
    match_value: __A | None = try_pick(chooser, xs)
    if match_value is None:
        return index_not_found()

    else:
        return value_1(match_value)


def reduce(folder: Callable[[_T, _T], _T], xs: IEnumerable_1[_T]) -> _T:
    with of_seq(xs) as e:

        def loop(acc_mut: _T | None = None, folder: Any = folder, xs: Any = xs) -> _T | None:
            while True:
                (acc,) = (acc_mut,)
                if e.System_Collections_IEnumerator_MoveNext():
                    acc_mut = folder(acc, e.System_Collections_Generic_IEnumerator_1_get_Current())
                    continue

                else:
                    return acc

                break

        if e.System_Collections_IEnumerator_MoveNext():
            return loop(e.System_Collections_Generic_IEnumerator_1_get_Current())

        else:
            raise Exception(SR_inputSequenceEmpty)


def reduce_back(folder: Callable[[_T, _T], _T], xs: IEnumerable_1[_T]) -> _T:
    arr: Array[_T] = to_array(xs)
    if len(arr) > 0:
        return reduce_back_1(folder, arr)

    else:
        raise Exception(SR_inputSequenceEmpty)


def replicate(n: int, x: __A) -> IEnumerable_1[__A]:
    def _arrow182(_arg: int, n: Any = n, x: Any = x) -> __A:
        return x

    return initialize(n, _arrow182)


def reverse(xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow183(__unit: None = None, xs: Any = xs) -> IEnumerable_1[_T]:
        return of_array(reverse_1(to_array(xs)))

    return delay(_arrow183)


def scan(folder: Callable[[_STATE, _T], _STATE], state: _STATE, xs: IEnumerable_1[_T]) -> IEnumerable_1[_STATE]:
    def _arrow184(__unit: None = None, folder: Any = folder, state: Any = state, xs: Any = xs) -> IEnumerable_1[_STATE]:
        acc: _STATE = state

        def mapping(x: _T | None = None) -> _STATE:
            nonlocal acc
            acc = folder(acc, x)
            return acc

        return concat(to_enumerable([singleton(state), map(mapping, xs)]))

    return delay(_arrow184)


def scan_back(folder: Callable[[_T, _STATE], _STATE], xs: IEnumerable_1[_T], state: _STATE) -> IEnumerable_1[_STATE]:
    def _arrow185(__unit: None = None, folder: Any = folder, xs: Any = xs, state: Any = state) -> IEnumerable_1[_STATE]:
        return of_array(scan_back_1(folder, to_array(xs), state, None))

    return delay(_arrow185)


def skip(count: int, source: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow186(__unit: None = None, count: Any = count, source: Any = source) -> IEnumerator[_T]:
        e: IEnumerator[_T] = of_seq(source)
        try:
            for _ in range(1, count + 1, 1):
                if not e.System_Collections_IEnumerator_MoveNext():
                    raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "source")

            def compensation(__unit: None = None) -> None:
                pass

            return Enumerator_enumerateThenFinally(compensation, e)

        except Exception as match_value:
            dispose_2(e)
            raise match_value

    return mk_seq(_arrow186)


def skip_while(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow187(__unit: None = None, predicate: Any = predicate, xs: Any = xs) -> IEnumerable_1[_T]:
        skipped: bool = True

        def f(x: _T | None = None) -> bool:
            nonlocal skipped
            if skipped:
                skipped = predicate(x)

            return not skipped

        return filter(f, xs)

    return delay(_arrow187)


def tail(xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    return skip(1, xs)


def take(count: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow188(__unit: None = None, count: Any = count, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow189(i: int, e: IEnumerator[_T], count: Any = count, xs: Any = xs) -> _T | None:
        if i < count:
            if e.System_Collections_IEnumerator_MoveNext():
                return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

            else:
                raise Exception((SR_notEnoughElements + "\\nParameter name: ") + "source")

        else:
            return None

    def _arrow190(e_1: IEnumerator[_T], count: Any = count, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow188, _arrow189, _arrow190)


def take_while(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow191(__unit: None = None, predicate: Any = predicate, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow192(e: IEnumerator[_T], predicate: Any = predicate, xs: Any = xs) -> _T | None:
        return (
            some(e.System_Collections_Generic_IEnumerator_1_get_Current())
            if (
                predicate(e.System_Collections_Generic_IEnumerator_1_get_Current())
                if e.System_Collections_IEnumerator_MoveNext()
                else False
            )
            else None
        )

    def _arrow193(e_1: IEnumerator[_T], predicate: Any = predicate, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate(_arrow191, _arrow192, _arrow193)


def truncate(count: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow194(__unit: None = None, count: Any = count, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow195(i: int, e: IEnumerator[_T], count: Any = count, xs: Any = xs) -> _T | None:
        return (
            some(e.System_Collections_Generic_IEnumerator_1_get_Current())
            if (e.System_Collections_IEnumerator_MoveNext() if (i < count) else False)
            else None
        )

    def _arrow196(e_1: IEnumerator[_T], count: Any = count, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow194, _arrow195, _arrow196)


def zip(xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2]) -> IEnumerable_1[tuple[_T1, _T2]]:
    def _arrow197(x: _T1, y: _T2, xs: Any = xs, ys: Any = ys) -> tuple[_T1, _T2]:
        return (x, y)

    return map2(_arrow197, xs, ys)


def zip3(xs: IEnumerable_1[_T1], ys: IEnumerable_1[_T2], zs: IEnumerable_1[_T3]) -> IEnumerable_1[tuple[_T1, _T2, _T3]]:
    def _arrow198(x: _T1, y: _T2, z: _T3, xs: Any = xs, ys: Any = ys, zs: Any = zs) -> tuple[_T1, _T2, _T3]:
        return (x, y, z)

    return map3(_arrow198, xs, ys, zs)


def collect(mapping: Callable[[_T], _COLLECTION], xs: IEnumerable_1[_T]) -> IEnumerable_1[Any]:
    def _arrow199(__unit: None = None, mapping: Any = mapping, xs: Any = xs) -> IEnumerable_1[_U]:
        return concat(map(mapping, xs))

    return delay(_arrow199)


def where(predicate: Callable[[_T], bool], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    return filter(predicate, xs)


def pairwise(xs: IEnumerable_1[_T]) -> IEnumerable_1[tuple[_T, _T]]:
    def _arrow200(__unit: None = None, xs: Any = xs) -> IEnumerable_1[tuple[_T, _T]]:
        return of_array(pairwise_1(to_array(xs)))

    return delay(_arrow200)


def split_into(chunks: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[Array[_T]]:
    def _arrow201(__unit: None = None, chunks: Any = chunks, xs: Any = xs) -> IEnumerable_1[Array[_T]]:
        return of_array(split_into_1(chunks, to_array(xs)))

    return delay(_arrow201)


def windowed(window_size: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[Array[_T]]:
    def _arrow202(__unit: None = None, window_size: Any = window_size, xs: Any = xs) -> IEnumerable_1[Array[_T]]:
        return of_array(windowed_1(window_size, to_array(xs)))

    return delay(_arrow202)


def transpose(xss: IEnumerable_1[Any]) -> IEnumerable_1[IEnumerable_1[Any]]:
    def _arrow203(__unit: None = None, xss: Any = xss) -> IEnumerable_1[IEnumerable_1[_T]]:
        return of_array(map_1(of_array, transpose_1(map_1(to_array, to_array(xss), None), None), None))

    return delay(_arrow203)


def sort_with(comparer: Callable[[_T, _T], int], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow204(__unit: None = None, comparer: Any = comparer, xs: Any = xs) -> IEnumerable_1[_T]:
        arr: Array[_T] = to_array(xs)
        sort_in_place_with(comparer, arr)
        return of_array(arr)

    return delay(_arrow204)


def sort(xs: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> IEnumerable_1[_T]:
    def _arrow205(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(x, y)

    return sort_with(_arrow205, xs)


def sort_by(projection: Callable[[_T], _U], xs: IEnumerable_1[_T], comparer: IComparer_1[_U]) -> IEnumerable_1[_T]:
    def _arrow206(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(projection(x), projection(y))

    return sort_with(_arrow206, xs)


def sort_descending(xs: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> IEnumerable_1[_T]:
    def _arrow207(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(x, y) * -1

    return sort_with(_arrow207, xs)


def sort_by_descending(
    projection: Callable[[_T], _U], xs: IEnumerable_1[_T], comparer: IComparer_1[_U]
) -> IEnumerable_1[_T]:
    def _arrow208(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> int:
        return comparer.Compare(projection(x), projection(y)) * -1

    return sort_with(_arrow208, xs)


def sum(xs: IEnumerable_1[_T], adder: IGenericAdder_1[_T]) -> _T:
    def _arrow209(acc: _T, x: _T, xs: Any = xs, adder: Any = adder) -> _T:
        return adder.Add(acc, x)

    return fold(_arrow209, adder.GetZero(), xs)


def sum_by(f: Callable[[_T], _U], xs: IEnumerable_1[_T], adder: IGenericAdder_1[_U]) -> _U:
    def _arrow210(acc: _U, x: _T, f: Any = f, xs: Any = xs, adder: Any = adder) -> _U:
        return adder.Add(acc, f(x))

    return fold(_arrow210, adder.GetZero(), xs)


def max_by(projection: Callable[[_T], _U], xs: IEnumerable_1[_T], comparer: IComparer_1[_U]) -> _T:
    def _arrow211(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> _T:
        return y if (comparer.Compare(projection(y), projection(x)) > 0) else x

    return reduce(_arrow211, xs)


def max(xs: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> _T:
    def _arrow212(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> _T:
        return y if (comparer.Compare(y, x) > 0) else x

    return reduce(_arrow212, xs)


def min_by(projection: Callable[[_T], _U], xs: IEnumerable_1[_T], comparer: IComparer_1[_U]) -> _T:
    def _arrow213(x: _T, y: _T, projection: Any = projection, xs: Any = xs, comparer: Any = comparer) -> _T:
        return x if (comparer.Compare(projection(y), projection(x)) > 0) else y

    return reduce(_arrow213, xs)


def min(xs: IEnumerable_1[_T], comparer: IComparer_1[_T]) -> _T:
    def _arrow214(x: _T, y: _T, xs: Any = xs, comparer: Any = comparer) -> _T:
        return x if (comparer.Compare(y, x) > 0) else y

    return reduce(_arrow214, xs)


def average(xs: IEnumerable_1[_T], averager: IGenericAverager_1[_T]) -> _T:
    count: int = 0

    def folder(acc: _T, x: _T, xs: Any = xs, averager: Any = averager) -> _T:
        nonlocal count
        count = (count + 1) or 0
        return averager.Add(acc, x)

    total: _T = fold(folder, averager.GetZero(), xs)
    if count == 0:
        raise Exception((SR_inputSequenceEmpty + "\\nParameter name: ") + "source")

    else:
        return averager.DivideByInt(total, count)


def average_by(f: Callable[[_T], _U], xs: IEnumerable_1[_T], averager: IGenericAverager_1[_U]) -> _U:
    count: int = 0

    def _arrow215(acc: _U, x: _T, f: Any = f, xs: Any = xs, averager: Any = averager) -> _U:
        nonlocal count
        count = (count + 1) or 0
        return averager.Add(acc, f(x))

    total: _U = fold(_arrow215, averager.GetZero(), xs)
    if count == 0:
        raise Exception((SR_inputSequenceEmpty + "\\nParameter name: ") + "source")

    else:
        return averager.DivideByInt(total, count)


def permute(f: Callable[[int], int], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    def _arrow216(__unit: None = None, f: Any = f, xs: Any = xs) -> IEnumerable_1[_T]:
        return of_array(permute_1(f, to_array(xs)))

    return delay(_arrow216)


def chunk_by_size(chunk_size: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[Array[_T]]:
    def _arrow217(__unit: None = None, chunk_size: Any = chunk_size, xs: Any = xs) -> IEnumerable_1[Array[_T]]:
        return of_array(chunk_by_size_1(chunk_size, to_array(xs)))

    return delay(_arrow217)


def insert_at(index: int, y: _T, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    is_done: bool = False
    if index < 0:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    def _arrow218(__unit: None = None, index: Any = index, y: Any = y, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow219(i: int, e: IEnumerator[_T], index: Any = index, y: Any = y, xs: Any = xs) -> _T | None:
        nonlocal is_done
        if e.System_Collections_IEnumerator_MoveNext() if (True if is_done else (i < index)) else False:
            return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

        elif i == index:
            is_done = True
            return some(y)

        else:
            if not is_done:
                raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

            return None

    def _arrow220(e_1: IEnumerator[_T], index: Any = index, y: Any = y, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow218, _arrow219, _arrow220)


def insert_many_at(index: int, ys: IEnumerable_1[_T], xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    status: int = -1
    if index < 0:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    def _arrow221(
        __unit: None = None, index: Any = index, ys: Any = ys, xs: Any = xs
    ) -> tuple[IEnumerator[_T], IEnumerator[_T]]:
        return (of_seq(xs), of_seq(ys))

    def _arrow222(
        i: int, tupled_arg: tuple[IEnumerator[_T], IEnumerator[_T]], index: Any = index, ys: Any = ys, xs: Any = xs
    ) -> _T | None:
        nonlocal status
        e1: IEnumerator[_T] = tupled_arg[0]
        e2: IEnumerator[_T] = tupled_arg[1]
        if i == index:
            status = 0

        inserted: _T | None
        if status == 0:
            if e2.System_Collections_IEnumerator_MoveNext():
                inserted = some(e2.System_Collections_Generic_IEnumerator_1_get_Current())

            else:
                status = 1
                inserted = None

        else:
            inserted = None

        if inserted is None:
            if e1.System_Collections_IEnumerator_MoveNext():
                return some(e1.System_Collections_Generic_IEnumerator_1_get_Current())

            else:
                if status < 1:
                    raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

                return None

        else:
            return some(value_1(inserted))

    def _arrow223(
        tupled_arg_1: tuple[IEnumerator[_T], IEnumerator[_T]], index: Any = index, ys: Any = ys, xs: Any = xs
    ) -> None:
        dispose_2(tupled_arg_1[0])
        dispose_2(tupled_arg_1[1])

    return generate_indexed(_arrow221, _arrow222, _arrow223)


def remove_at(index: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    is_done: bool = False
    if index < 0:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    def _arrow224(__unit: None = None, index: Any = index, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow225(i: int, e: IEnumerator[_T], index: Any = index, xs: Any = xs) -> _T | None:
        nonlocal is_done
        if e.System_Collections_IEnumerator_MoveNext() if (True if is_done else (i < index)) else False:
            return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

        elif e.System_Collections_IEnumerator_MoveNext() if (i == index) else False:
            is_done = True
            return (
                some(e.System_Collections_Generic_IEnumerator_1_get_Current())
                if e.System_Collections_IEnumerator_MoveNext()
                else None
            )

        else:
            if not is_done:
                raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

            return None

    def _arrow226(e_1: IEnumerator[_T], index: Any = index, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow224, _arrow225, _arrow226)


def remove_many_at(index: int, count: int, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    if index < 0:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    def _arrow227(__unit: None = None, index: Any = index, count: Any = count, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow228(i: int, e: IEnumerator[_T], index: Any = index, count: Any = count, xs: Any = xs) -> _T | None:
        if i < index:
            if e.System_Collections_IEnumerator_MoveNext():
                return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

            else:
                raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

        else:
            if i == index:
                for _ in range(1, count + 1, 1):
                    if not e.System_Collections_IEnumerator_MoveNext():
                        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "count")

            return (
                some(e.System_Collections_Generic_IEnumerator_1_get_Current())
                if e.System_Collections_IEnumerator_MoveNext()
                else None
            )

    def _arrow229(e_1: IEnumerator[_T], index: Any = index, count: Any = count, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow227, _arrow228, _arrow229)


def update_at(index: int, y: _T, xs: IEnumerable_1[_T]) -> IEnumerable_1[_T]:
    is_done: bool = False
    if index < 0:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    def _arrow230(__unit: None = None, index: Any = index, y: Any = y, xs: Any = xs) -> IEnumerator[_T]:
        return of_seq(xs)

    def _arrow231(i: int, e: IEnumerator[_T], index: Any = index, y: Any = y, xs: Any = xs) -> _T | None:
        nonlocal is_done
        if e.System_Collections_IEnumerator_MoveNext() if (True if is_done else (i < index)) else False:
            return some(e.System_Collections_Generic_IEnumerator_1_get_Current())

        elif e.System_Collections_IEnumerator_MoveNext() if (i == index) else False:
            is_done = True
            return some(y)

        else:
            if not is_done:
                raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

            return None

    def _arrow232(e_1: IEnumerator[_T], index: Any = index, y: Any = y, xs: Any = xs) -> None:
        dispose_2(e_1)

    return generate_indexed(_arrow230, _arrow231, _arrow232)


__all__ = [
    "SR_enumerationAlreadyFinished",
    "SR_enumerationNotStarted",
    "SR_inputSequenceEmpty",
    "SR_inputSequenceTooLong",
    "SR_keyNotFoundAlt",
    "SR_notEnoughElements",
    "SR_resetNotSupported",
    "Enumerator_noReset",
    "Enumerator_notStarted",
    "Enumerator_alreadyFinished",
    "Enumerator_Seq_reflection",
    "Enumerator_FromFunctions_1_reflection",
    "Enumerator_cast",
    "Enumerator_concat",
    "Enumerator_enumerateThenFinally",
    "Enumerator_generateWhileSome",
    "Enumerator_unfold",
    "index_not_found",
    "check_non_null",
    "mk_seq",
    "of_seq",
    "delay",
    "concat",
    "unfold",
    "empty",
    "singleton",
    "of_array",
    "to_array",
    "of_list",
    "to_list",
    "generate",
    "generate_indexed",
    "append",
    "cast",
    "choose",
    "compare_with",
    "contains",
    "enumerate_from_functions",
    "enumerate_then_finally",
    "enumerate_using",
    "enumerate_while",
    "filter",
    "exists",
    "exists2",
    "exactly_one",
    "try_exactly_one",
    "try_find",
    "find",
    "try_find_back",
    "find_back",
    "try_find_index",
    "find_index",
    "try_find_index_back",
    "find_index_back",
    "fold",
    "fold_back",
    "fold2",
    "fold_back2",
    "for_all",
    "for_all2",
    "try_head",
    "head",
    "initialize",
    "initialize_infinite",
    "is_empty",
    "try_item",
    "item",
    "iterate",
    "iterate2",
    "iterate_indexed",
    "iterate_indexed2",
    "try_last",
    "last",
    "length",
    "map",
    "map_indexed",
    "indexed",
    "map2",
    "map_indexed2",
    "map3",
    "read_only",
    "CachedSeq_1_reflection",
    "CachedSeq_1__Clear",
    "cache",
    "all_pairs",
    "map_fold",
    "map_fold_back",
    "try_pick",
    "pick",
    "reduce",
    "reduce_back",
    "replicate",
    "reverse",
    "scan",
    "scan_back",
    "skip",
    "skip_while",
    "tail",
    "take",
    "take_while",
    "truncate",
    "zip",
    "zip3",
    "collect",
    "where",
    "pairwise",
    "split_into",
    "windowed",
    "transpose",
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
    "insert_at",
    "insert_many_at",
    "remove_at",
    "remove_many_at",
    "update_at",
]
