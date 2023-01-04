from __future__ import annotations
from typing import (TypeVar, Any, Optional, Tuple, Callable, Generic)
from .map_util import (try_get_value, get_item_from_dict)
from .reflection import (TypeInfo, class_type)
from .resize_array import find_index
from .seq import (concat, iterate_indexed, to_array, delay, map, iterate)
from .string import format
from .types import (Array, FSharpRef)
from .util import (get_enumerator, IEnumerator, to_iterator, equals, ignore, IEnumerable_1, ICollection, IEqualityComparer_1, dispose)

_KEY = TypeVar("_KEY")

_VALUE = TypeVar("_VALUE")

def _expr32(gen0: TypeInfo, gen1: TypeInfo) -> TypeInfo:
    return class_type("Fable.Collections.Dictionary", [gen0, gen1], Dictionary)


class Dictionary(Generic[_KEY, _VALUE]):
    def __init__(self, pairs: IEnumerable_1[Any], comparer: IEqualityComparer_1[Any]) -> None:
        this: FSharpRef[Dictionary[_KEY, _VALUE]] = FSharpRef(None)
        self.comparer: IEqualityComparer_1[Any] = comparer
        this.contents = self
        self.hash_map: Any = dict([])
        self.init_00409: int = 1
        with get_enumerator(pairs) as enumerator:
            while enumerator.System_Collections_IEnumerator_MoveNext():
                pair: Any = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
                Dictionary__Add_5BDDA1(this.contents, pair[0], pair[1])

    @property
    def Symbol_toStringTag(self, __unit: None=None) -> str:
        return "Dictionary"

    def to_json(self, __unit: None=None) -> Any:
        this: Dictionary[_KEY, _VALUE] = self
        return list(this)

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        this: Dictionary[_KEY, _VALUE] = self
        return get_enumerator(this)

    def GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        this: Dictionary[_KEY, _VALUE] = self
        return get_enumerator(concat(this.hash_map.values()))

    def __iter__(self) -> IEnumerator[Any]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_Generic_ICollection_1_Add2B595(self, item: Any) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__Add_5BDDA1(this, item[0], item[1])

    def System_Collections_Generic_ICollection_1_Clear(self, __unit: None=None) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__Clear(this)

    def System_Collections_Generic_ICollection_1_Contains2B595(self, item: Any) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        match_value: Optional[Any] = Dictionary__TryFind_2B595(this, item[0])
        (pattern_matching_result,) = (None,)
        if match_value is not None:
            if equals(match_value[1], item[1]):
                pattern_matching_result = 0

            else: 
                pattern_matching_result = 1


        else: 
            pattern_matching_result = 1

        if pattern_matching_result == 0:
            return True

        elif pattern_matching_result == 1:
            return False


    def System_Collections_Generic_ICollection_1_CopyToZ3B4C077E(self, array: Array[Any], array_index: int) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        def action(i: int, e: Any) -> None:
            array[array_index + i] = e

        iterate_indexed(action, this)

    def __len__(self, __unit: None=None) -> int:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__get_Count(this)

    def System_Collections_Generic_ICollection_1_get_IsReadOnly(self, __unit: None=None) -> bool:
        return False

    def System_Collections_Generic_ICollection_1_Remove2B595(self, item: Any) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        match_value: Optional[Any] = Dictionary__TryFind_2B595(this, item[0])
        if match_value is not None:
            if equals(match_value[1], item[1]):
                ignore(Dictionary__Remove_2B595(this, item[0]))

            return True

        else: 
            return False


    def System_Collections_Generic_IDictionary_2_Add5BDDA1(self, key: _KEY, value: _VALUE) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__Add_5BDDA1(this, key, value)

    def System_Collections_Generic_IDictionary_2_ContainsKey2B595(self, key: Optional[_KEY]=None) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__ContainsKey_2B595(this, key)

    def System_Collections_Generic_IDictionary_2_get_Item2B595(self, key: Optional[_KEY]=None) -> _VALUE:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__get_Item_2B595(this, key)

    def System_Collections_Generic_IDictionary_2_set_Item5BDDA1(self, key: _KEY, v: _VALUE) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__set_Item_5BDDA1(this, key, v)

    def System_Collections_Generic_IDictionary_2_get_Keys(self, __unit: None=None) -> ICollection[_KEY]:
        this: Dictionary[_KEY, _VALUE] = self
        def _arrow29(__unit: None=None) -> IEnumerable_1[_KEY]:
            def _arrow28(pair: Any) -> _KEY:
                return pair[0]

            return map(_arrow28, this)

        return to_array(delay(_arrow29))

    def System_Collections_Generic_IDictionary_2_Remove2B595(self, key: Optional[_KEY]=None) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__Remove_2B595(this, key)

    def System_Collections_Generic_IDictionary_2_TryGetValue6DC89625(self, key: _KEY, value: FSharpRef[_VALUE]) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        match_value: Optional[Any] = Dictionary__TryFind_2B595(this, key)
        if match_value is not None:
            pair: Any = match_value
            value.contents = pair[1]
            return True

        else: 
            return False


    def System_Collections_Generic_IDictionary_2_get_Values(self, __unit: None=None) -> ICollection[_VALUE]:
        this: Dictionary[_KEY, _VALUE] = self
        def _arrow31(__unit: None=None) -> IEnumerable_1[_VALUE]:
            def _arrow30(pair: Any) -> _VALUE:
                return pair[1]

            return map(_arrow30, this)

        return to_array(delay(_arrow31))

    @property
    def size(self, __unit: None=None) -> int:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__get_Count(this)

    def clear(self, __unit: None=None) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__Clear(this)

    def delete(self, k: Optional[_KEY]=None) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__Remove_2B595(this, k)

    def entries(self, __unit: None=None) -> IEnumerable_1[Tuple[_KEY, _VALUE]]:
        this: Dictionary[_KEY, _VALUE] = self
        def mapping(p: Any) -> Tuple[_KEY, _VALUE]:
            return (p[0], p[1])

        return map(mapping, this)

    def __getitem__(self, k: Optional[_KEY]=None) -> _VALUE:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__get_Item_2B595(this, k)

    def __contains__(self, k: Optional[_KEY]=None) -> bool:
        this: Dictionary[_KEY, _VALUE] = self
        return Dictionary__ContainsKey_2B595(this, k)

    def keys(self, __unit: None=None) -> IEnumerable_1[_KEY]:
        this: Dictionary[_KEY, _VALUE] = self
        def mapping(p: Any) -> _KEY:
            return p[0]

        return map(mapping, this)

    def __setitem__(self, k: _KEY, v: _VALUE) -> Map_2[_KEY, _VALUE]:
        this: Dictionary[_KEY, _VALUE] = self
        Dictionary__set_Item_5BDDA1(this, k, v)
        return this

    def values(self, __unit: None=None) -> IEnumerable_1[_VALUE]:
        this: Dictionary[_KEY, _VALUE] = self
        def mapping(p: Any) -> _VALUE:
            return p[1]

        return map(mapping, this)

    def for_each(self, f: Callable[[_VALUE, _KEY, Map_2[_KEY, _VALUE]], None], this_arg: Optional[Any]=None) -> None:
        this: Dictionary[_KEY, _VALUE] = self
        def action(p: Any) -> None:
            f(p[1], p[0], this)

        iterate(action, this)


Dictionary_reflection = _expr32

def Dictionary__ctor_6623D9B3(pairs: IEnumerable_1[Any], comparer: IEqualityComparer_1[Any]) -> Dictionary[_KEY, _VALUE]:
    return Dictionary(pairs, comparer)


def Dictionary__TryFindIndex_2B595(this: Dictionary[_KEY, Any], k: _KEY) -> Tuple[bool, int, int]:
    h: int = this.comparer.GetHashCode(k) or 0
    match_value: Tuple[bool, Array[Any]]
    out_arg: Array[Any] = None
    def _arrow33(__unit: None=None, this: Any=this, k: _KEY=k) -> Array[Any]:
        return out_arg

    def _arrow34(v: Array[Any], this: Any=this, k: _KEY=k) -> None:
        nonlocal out_arg
        out_arg = v

    match_value = (try_get_value(this.hash_map, h, FSharpRef(_arrow33, _arrow34)), out_arg)
    if match_value[0]:
        def _arrow35(pair: Any, this: Any=this, k: _KEY=k) -> bool:
            return this.comparer.Equals(k, pair[0])

        return (True, h, find_index(_arrow35, match_value[1]))

    else: 
        return (False, h, -1)



def Dictionary__TryFind_2B595(this: Dictionary[_KEY, _VALUE], k: _KEY) -> Optional[Any]:
    match_value: Tuple[bool, int, int] = Dictionary__TryFindIndex_2B595(this, k)
    (pattern_matching_result,) = (None,)
    if match_value[0]:
        if match_value[2] > -1:
            pattern_matching_result = 0

        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        return get_item_from_dict(this.hash_map, match_value[1])[match_value[2]]

    elif pattern_matching_result == 1:
        return None



def Dictionary__get_Comparer(this: Dictionary[_KEY, Any]) -> IEqualityComparer_1[Any]:
    return this.comparer


def Dictionary__Clear(this: Dictionary[Any, Any]) -> None:
    this.hash_map.clear()


def Dictionary__get_Count(this: Dictionary[Any, Any]) -> int:
    count: int = 0
    enumerator: Any = get_enumerator(this.hash_map.values())
    try: 
        while enumerator.System_Collections_IEnumerator_MoveNext():
            pairs: Array[Any] = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
            count = (count + len(pairs)) or 0

    finally: 
        dispose(enumerator)

    return count


def Dictionary__get_Item_2B595(this: Dictionary[_KEY, _VALUE], k: _KEY) -> _VALUE:
    match_value: Optional[Any] = Dictionary__TryFind_2B595(this, k)
    if match_value is not None:
        return match_value[1]

    else: 
        raise Exception("The item was not found in collection")



def Dictionary__set_Item_5BDDA1(this: Dictionary[_KEY, _VALUE], k: _KEY, v: _VALUE) -> None:
    match_value: Tuple[bool, int, int] = Dictionary__TryFindIndex_2B595(this, k)
    if match_value[0]:
        if match_value[2] > -1:
            get_item_from_dict(this.hash_map, match_value[1])[match_value[2]] = (k, v)

        else: 
            value: None = (get_item_from_dict(this.hash_map, match_value[1]).append((k, v)))
            ignore()


    else: 
        this.hash_map[match_value[1]] = [(k, v)]



def Dictionary__Add_5BDDA1(this: Dictionary[_KEY, _VALUE], k: _KEY, v: _VALUE) -> None:
    match_value: Tuple[bool, int, int] = Dictionary__TryFindIndex_2B595(this, k)
    if match_value[0]:
        if match_value[2] > -1:
            raise Exception(format("An item with the same key has already been added. Key: {0}", k))

        else: 
            value: None = (get_item_from_dict(this.hash_map, match_value[1]).append((k, v)))
            ignore()


    else: 
        this.hash_map[match_value[1]] = [(k, v)]



def Dictionary__ContainsKey_2B595(this: Dictionary[_KEY, Any], k: _KEY) -> bool:
    match_value: Tuple[bool, int, int] = Dictionary__TryFindIndex_2B595(this, k)
    (pattern_matching_result,) = (None,)
    if match_value[0]:
        if match_value[2] > -1:
            pattern_matching_result = 0

        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        return True

    elif pattern_matching_result == 1:
        return False



def Dictionary__Remove_2B595(this: Dictionary[_KEY, Any], k: _KEY) -> bool:
    match_value: Tuple[bool, int, int] = Dictionary__TryFindIndex_2B595(this, k)
    (pattern_matching_result,) = (None,)
    if match_value[0]:
        if match_value[2] > -1:
            pattern_matching_result = 0

        else: 
            pattern_matching_result = 1


    else: 
        pattern_matching_result = 1

    if pattern_matching_result == 0:
        get_item_from_dict(this.hash_map, match_value[1]).pop(match_value[2])
        return True

    elif pattern_matching_result == 1:
        return False



__all__ = ["Dictionary_reflection", "Dictionary__TryFindIndex_2B595", "Dictionary__TryFind_2B595", "Dictionary__get_Comparer", "Dictionary__Clear", "Dictionary__get_Count", "Dictionary__get_Item_2B595", "Dictionary__set_Item_5BDDA1", "Dictionary__Add_5BDDA1", "Dictionary__ContainsKey_2B595", "Dictionary__Remove_2B595"]

