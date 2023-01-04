from typing import (TypeVar, Any, Optional, Callable, Tuple)
from .list import FSharpList
from .map_util import (add_to_set, try_get_value, get_item_from_dict, add_to_dict)
from .mutable_map import Dictionary
from .mutable_set import HashSet
from .seq import (delay, filter, map, to_array, to_list)
from .types import (Array, FSharpRef)
from .util import (IEnumerable_1, IEqualityComparer_1, get_enumerator)

_T = TypeVar("_T")

_KEY = TypeVar("_KEY")

def distinct(xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[_T]:
    def _arrow114(__unit: None=None, xs: IEnumerable_1[Any]=xs, comparer: Any=comparer) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet([], comparer)
        def predicate(x: Optional[_T]=None) -> bool:
            return add_to_set(x, hash_set)

        return filter(predicate, xs)

    return delay(_arrow114)


def distinct_by(projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[_T]:
    def _arrow115(__unit: None=None, projection: Any=projection, xs: IEnumerable_1[Any]=xs, comparer: Any=comparer) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet([], comparer)
        def predicate(x: Optional[_T]=None) -> bool:
            return add_to_set(projection(x), hash_set)

        return filter(predicate, xs)

    return delay(_arrow115)


def except_(items_to_exclude: IEnumerable_1[_T], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[_T]:
    def _arrow116(__unit: None=None, items_to_exclude: IEnumerable_1[Any]=items_to_exclude, xs: IEnumerable_1[Any]=xs, comparer: Any=comparer) -> IEnumerable_1[_T]:
        hash_set: Any = HashSet(items_to_exclude, comparer)
        def predicate(x: Optional[_T]=None) -> bool:
            return add_to_set(x, hash_set)

        return filter(predicate, xs)

    return delay(_arrow116)


def count_by(projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[Tuple[_KEY, int]]:
    def _arrow120(__unit: None=None, projection: Any=projection, xs: IEnumerable_1[Any]=xs, comparer: Any=comparer) -> IEnumerable_1[Tuple[_KEY, int]]:
        dict_1: Any = Dictionary([], comparer)
        keys: Array[_KEY] = []
        with get_enumerator(xs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                key: _KEY = projection(enumerator.System_Collections_Generic_IEnumerator_1_get_Current())
                match_value: Tuple[bool, int]
                out_arg: int = 0
                def _arrow117(__unit: None=None) -> int:
                    return out_arg

                def _arrow118(v: int) -> None:
                    nonlocal out_arg
                    out_arg = v or 0

                match_value = (try_get_value(dict_1, key, FSharpRef(_arrow117, _arrow118)), out_arg)
                if match_value[0]:
                    dict_1[key] = match_value[1] + 1

                else: 
                    dict_1[key] = 1
                    (keys.append(key))

        def _arrow119(key_1: Optional[_KEY]=None) -> Tuple[_KEY, int]:
            return (key_1, get_item_from_dict(dict_1, key_1))

        return map(_arrow119, keys)

    return delay(_arrow120)


def group_by(projection: Callable[[_T], _KEY], xs: IEnumerable_1[_T], comparer: IEqualityComparer_1[Any]) -> IEnumerable_1[Tuple[_KEY, IEnumerable_1[_T]]]:
    def _arrow124(__unit: None=None, projection: Any=projection, xs: IEnumerable_1[Any]=xs, comparer: Any=comparer) -> IEnumerable_1[Tuple[_KEY, IEnumerable_1[_T]]]:
        dict_1: Any = Dictionary([], comparer)
        keys: Array[_KEY] = []
        with get_enumerator(xs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                x: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
                key: _KEY = projection(x)
                match_value: Tuple[bool, Array[_T]]
                out_arg: Array[_T] = None
                def _arrow121(__unit: None=None) -> Array[_T]:
                    return out_arg

                def _arrow122(v: Array[_T]) -> None:
                    nonlocal out_arg
                    out_arg = v

                match_value = (try_get_value(dict_1, key, FSharpRef(_arrow121, _arrow122)), out_arg)
                if match_value[0]:
                    (match_value[1].append(x))

                else: 
                    add_to_dict(dict_1, key, [x])
                    (keys.append(key))

        def _arrow123(key_1: Optional[_KEY]=None) -> Tuple[_KEY, IEnumerable_1[_T]]:
            return (key_1, get_item_from_dict(dict_1, key_1))

        return map(_arrow123, keys)

    return delay(_arrow124)


def Array_distinct(xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(distinct(xs, comparer))


def Array_distinctBy(projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(distinct_by(projection, xs, comparer))


def Array_except(items_to_exclude: IEnumerable_1[_T], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[_T]:
    return to_array(except_(items_to_exclude, xs, comparer))


def Array_countBy(projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[Tuple[_KEY, int]]:
    return to_array(count_by(projection, xs, comparer))


def Array_groupBy(projection: Callable[[_T], _KEY], xs: Array[_T], comparer: IEqualityComparer_1[Any]) -> Array[Tuple[_KEY, Array[_T]]]:
    def mapping(tupled_arg: Tuple[_KEY, IEnumerable_1[_T]], projection: Any=projection, xs: Array[Any]=xs, comparer: Any=comparer) -> Tuple[_KEY, Array[_T]]:
        return (tupled_arg[0], to_array(tupled_arg[1]))

    return to_array(map(mapping, group_by(projection, xs, comparer)))


def List_distinct(xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[_T]:
    return to_list(distinct(xs, comparer))


def List_distinctBy(projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[_T]:
    return to_list(distinct_by(projection, xs, comparer))


def List_except(items_to_exclude: IEnumerable_1[_T], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[_T]:
    return to_list(except_(items_to_exclude, xs, comparer))


def List_countBy(projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[Tuple[_KEY, int]]:
    return to_list(count_by(projection, xs, comparer))


def List_groupBy(projection: Callable[[_T], _KEY], xs: FSharpList[_T], comparer: IEqualityComparer_1[Any]) -> FSharpList[Tuple[_KEY, FSharpList[_T]]]:
    def mapping(tupled_arg: Tuple[_KEY, IEnumerable_1[_T]], projection: Any=projection, xs: FSharpList[Any]=xs, comparer: Any=comparer) -> Tuple[_KEY, FSharpList[_T]]:
        return (tupled_arg[0], to_list(tupled_arg[1]))

    return to_list(map(mapping, group_by(projection, xs, comparer)))


__all__ = ["distinct", "distinct_by", "except_", "count_by", "group_by", "Array_distinct", "Array_distinctBy", "Array_except", "Array_countBy", "Array_groupBy", "List_distinct", "List_distinctBy", "List_except", "List_countBy", "List_groupBy"]

