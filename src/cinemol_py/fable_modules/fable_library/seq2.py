from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

from .list import FSharpList
from .map_util import add_to_dict, add_to_set, get_item_from_dict, try_get_value
from .mutable_map import Dictionary
from .mutable_set import HashSet
from .seq import delay, filter, map, to_array, to_list
from .types import Array, FSharpRef
from .util import IEnumerable_1, IEqualityComparer_1, get_enumerator


_T = TypeVar("_T")

_KEY = TypeVar("_KEY")


def distinct(xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[_T]:
    def _arrow117(__unit: None = None, xs: Any = xs, comparer: Any = comparer) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet([], comparer)

        def predicate(x: _T | None = None) -> bool:
            return add_to_set(x, hash_set)

        return filter(predicate, xs)

    return delay(_arrow117)


def distinct_by(
    projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]
) -> IEnumerable_1[_T]:
    def _arrow118(
        __unit: None = None, projection: Any = projection, xs: Any = xs, comparer: Any = comparer
    ) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet([], comparer)

        def predicate(x: _T | None = None) -> bool:
            return add_to_set(projection(x), hash_set)

        return filter(predicate, xs)

    return delay(_arrow118)


def except_(
    items_to_exclude: IEnumerable_1[_T], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]
) -> IEnumerable_1[_T]:
    def _arrow119(
        __unit: None = None, items_to_exclude: Any = items_to_exclude, xs: Any = xs, comparer: Any = comparer
    ) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet(items_to_exclude, comparer)

        def predicate(x: _T | None = None) -> bool:
            return add_to_set(x, hash_set)

        return filter(predicate, xs)

    return delay(_arrow119)


def count_by(
    projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]
) -> IEnumerable_1[tuple[_KEY, int]]:
    def _arrow123(
        __unit: None = None, projection: Any = projection, xs: Any = xs, comparer: Any = comparer
    ) -> IEnumerable_1[tuple[_KEY, int]]:
        dict_1: Any = Dictionary([], comparer)
        keys: Array[_KEY] = []
        with get_enumerator(xs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                key: _KEY = projection(enumerator.System_Collections_Generic_IEnumerator_1_get_Current())
                match_value: tuple[bool, int]
                out_arg: int = 0

                def _arrow120(__unit: None = None) -> int:
                    return out_arg

                def _arrow121(v: int) -> None:
                    nonlocal out_arg
                    out_arg = v or 0

                match_value = (try_get_value(dict_1, key, FSharpRef(_arrow120, _arrow121)), out_arg)
                if match_value[0]:
                    dict_1[key] = match_value[1] + 1

                else:
                    dict_1[key] = 1
                    (keys.append(key))

        def _arrow122(key_1: _KEY | None = None) -> tuple[_KEY, int]:
            return (key_1, get_item_from_dict(dict_1, key_1))

        return map(_arrow122, keys)

    return delay(_arrow123)


def group_by(
    projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]
) -> IEnumerable_1[tuple[_KEY, IEnumerable_1[_T]]]:
    def _arrow127(
        __unit: None = None, projection: Any = projection, xs: Any = xs, comparer: Any = comparer
    ) -> IEnumerable_1[tuple[_KEY, IEnumerable_1[_T]]]:
        dict_1: Any = Dictionary([], comparer)
        keys: Array[_KEY] = []
        with get_enumerator(xs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                x: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
                key: _KEY = projection(x)
                match_value: tuple[bool, Array[_T]]
                out_arg: Array[_T] = None

                def _arrow124(__unit: None = None) -> Array[_T]:
                    return out_arg

                def _arrow125(v: Array[_T]) -> None:
                    nonlocal out_arg
                    out_arg = v

                match_value = (try_get_value(dict_1, key, FSharpRef(_arrow124, _arrow125)), out_arg)
                if match_value[0]:
                    (match_value[1].append(x))

                else:
                    add_to_dict(dict_1, key, [x])
                    (keys.append(key))

        def _arrow126(key_1: _KEY | None = None) -> tuple[_KEY, IEnumerable_1[_T]]:
            return (key_1, get_item_from_dict(dict_1, key_1))

        return map(_arrow126, keys)

    return delay(_arrow127)


def Array_distinct(xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(distinct(xs, comparer))


def Array_distinctBy(projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(distinct_by(projection, xs, comparer))


def Array_except(items_to_exclude: IEnumerable_1[_T], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(except_(items_to_exclude, xs, comparer))


def Array_countBy(
    projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]
) -> Array[tuple[_KEY, int]]:
    return to_array(count_by(projection, xs, comparer))


def Array_groupBy(
    projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]
) -> Array[tuple[_KEY, Array[_T]]]:
    def mapping(
        tupled_arg: tuple[_KEY, IEnumerable_1[_T]], projection: Any = projection, xs: Any = xs, comparer: Any = comparer
    ) -> tuple[_KEY, Array[_T]]:
        return (tupled_arg[0], to_array(tupled_arg[1]))

    return to_array(map(mapping, group_by(projection, xs, comparer)))


def List_distinct(xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[_T]:
    return to_list(distinct(xs, comparer))


def List_distinctBy(
    projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]
) -> FSharpList[_T]:
    return to_list(distinct_by(projection, xs, comparer))


def List_except(
    items_to_exclude: IEnumerable_1[_T], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]
) -> FSharpList[_T]:
    return to_list(except_(items_to_exclude, xs, comparer))


def List_countBy(
    projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]
) -> FSharpList[tuple[_KEY, int]]:
    return to_list(count_by(projection, xs, comparer))


def List_groupBy(
    projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]
) -> FSharpList[tuple[_KEY, FSharpList[_T]]]:
    def mapping(
        tupled_arg: tuple[_KEY, IEnumerable_1[_T]], projection: Any = projection, xs: Any = xs, comparer: Any = comparer
    ) -> tuple[_KEY, FSharpList[_T]]:
        return (tupled_arg[0], to_list(tupled_arg[1]))

    return to_list(map(mapping, group_by(projection, xs, comparer)))


__all__ = [
    "distinct",
    "distinct_by",
    "except_",
    "count_by",
    "group_by",
    "Array_distinct",
    "Array_distinctBy",
    "Array_except",
    "Array_countBy",
    "Array_groupBy",
    "List_distinct",
    "List_distinctBy",
    "List_except",
    "List_countBy",
    "List_groupBy",
]
