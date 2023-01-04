from __future__ import annotations
import builtins
import functools
from math import ceil
from typing import (Any, TypeVar, Callable, Optional, Tuple)
from .double import divide
from .option import (some, default_arg, map as map_1, value as value_1)
from .types import Array
from .util import (max as max_1, compare_primitives, IEnumerable_1, IEqualityComparer_1, copy_to_array, ignore, get_enumerator, equals as equals_1, IComparer_1, min as min_1)
from .native import (Cons_1, Helpers_allocateArrayFromCons, Helpers_fillImpl, Helpers_spliceImpl)
from .global_ import (SR_indexOutOfBounds, IGenericAdder_1, IGenericAverager_1)

_T = TypeVar("_T")

_U = TypeVar("_U")

_T1 = TypeVar("_T1")

_T2 = TypeVar("_T2")

_T3 = TypeVar("_T3")

_STATE = TypeVar("_STATE")

_RESULT = TypeVar("_RESULT")

__A = TypeVar("__A")

__B = TypeVar("__B")

_A = TypeVar("_A")

_B = TypeVar("_B")

__C = TypeVar("__C")

def index_not_found(__unit: None=None) -> Any:
    raise Exception("An index satisfying the predicate was not found in the collection.")


def different_lengths(__unit: None=None) -> Any:
    raise Exception("Arrays had different lengths")


def append(array1: Array[_T], array2: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    len1: int = len(array1) or 0
    len2: int = len(array2) or 0
    new_array: Array[_T] = Helpers_allocateArrayFromCons(cons, len1 + len2)
    for i in range(0, (len1 - 1) + 1, 1):
        new_array[i] = array1[i]
    for i_1 in range(0, (len2 - 1) + 1, 1):
        new_array[i_1 + len1] = array2[i_1]
    return new_array


def filter(predicate: Callable[[_T], bool], array: Array[_T]) -> Array[_T]:
    return list(builtins.filter(predicate, array))


def fill(target: Array[_T], target_index: int, count: int, value: _T) -> Array[_T]:
    return Helpers_fillImpl(target, value, target_index, count)


def get_sub_array(array: Array[_T], start: int, count: int) -> Array[_T]:
    return array[start:start+count]


def last(array: Array[_T]) -> _T:
    if len(array) == 0:
        raise Exception("The input array was empty\\nParameter name: array")

    return array[len(array) - 1]


def try_last(array: Array[_T]) -> Optional[_T]:
    if len(array) == 0:
        return None

    else: 
        return some(array[len(array) - 1])



def map_indexed(f: Callable[[int, _T], _U], source: Array[_T], cons: Cons_1[_U]) -> Array[_U]:
    len_1: int = len(source) or 0
    target: Array[_U] = Helpers_allocateArrayFromCons(cons, len_1)
    for i in range(0, (len_1 - 1) + 1, 1):
        target[i] = f(i, source[i])
    return target


def map(f: Callable[[_T], _U], source: Array[_T], cons: Cons_1[_U]) -> Array[_U]:
    len_1: int = len(source) or 0
    target: Array[_U] = Helpers_allocateArrayFromCons(cons, len_1)
    for i in range(0, (len_1 - 1) + 1, 1):
        target[i] = f(source[i])
    return target


def map_indexed2(f: Callable[[int, _T1, _T2], _U], source1: Array[_T1], source2: Array[_T2], cons: Cons_1[_U]) -> Array[_U]:
    if len(source1) != len(source2):
        raise Exception("Arrays had different lengths")

    result: Array[_U] = Helpers_allocateArrayFromCons(cons, len(source1))
    for i in range(0, (len(source1) - 1) + 1, 1):
        result[i] = f(i, source1[i], source2[i])
    return result


def map2(f: Callable[[_T1, _T2], _U], source1: Array[_T1], source2: Array[_T2], cons: Cons_1[_U]) -> Array[_U]:
    if len(source1) != len(source2):
        raise Exception("Arrays had different lengths")

    result: Array[_U] = Helpers_allocateArrayFromCons(cons, len(source1))
    for i in range(0, (len(source1) - 1) + 1, 1):
        result[i] = f(source1[i], source2[i])
    return result


def map_indexed3(f: Callable[[int, _T1, _T2, _T3], _U], source1: Array[_T1], source2: Array[_T2], source3: Array[_T3], cons: Cons_1[_U]) -> Array[_U]:
    if True if (len(source1) != len(source2)) else (len(source2) != len(source3)):
        raise Exception("Arrays had different lengths")

    result: Array[_U] = Helpers_allocateArrayFromCons(cons, len(source1))
    for i in range(0, (len(source1) - 1) + 1, 1):
        result[i] = f(i, source1[i], source2[i], source3[i])
    return result


def map3(f: Callable[[_T1, _T2, _T3], _U], source1: Array[_T1], source2: Array[_T2], source3: Array[_T3], cons: Cons_1[_U]) -> Array[_U]:
    if True if (len(source1) != len(source2)) else (len(source2) != len(source3)):
        raise Exception("Arrays had different lengths")

    result: Array[_U] = Helpers_allocateArrayFromCons(cons, len(source1))
    for i in range(0, (len(source1) - 1) + 1, 1):
        result[i] = f(source1[i], source2[i], source3[i])
    return result


def map_fold(mapping: Callable[[_STATE, _T], Tuple[_RESULT, _STATE]], state: _STATE, array: Array[_T], cons: Cons_1[_RESULT]) -> Tuple[Array[_RESULT], _STATE]:
    match_value: int = len(array) or 0
    if match_value == 0:
        return ([], state)

    else: 
        acc: _STATE = state
        res: Array[_RESULT] = Helpers_allocateArrayFromCons(cons, match_value)
        for i in range(0, (len(array) - 1) + 1, 1):
            pattern_input: Tuple[_RESULT, _STATE] = mapping(acc, array[i])
            res[i] = pattern_input[0]
            acc = pattern_input[1]
        return (res, acc)



def map_fold_back(mapping: Callable[[_T, _STATE], Tuple[_RESULT, _STATE]], array: Array[_T], state: _STATE, cons: Cons_1[_RESULT]) -> Tuple[Array[_RESULT], _STATE]:
    match_value: int = len(array) or 0
    if match_value == 0:
        return ([], state)

    else: 
        acc: _STATE = state
        res: Array[_RESULT] = Helpers_allocateArrayFromCons(cons, match_value)
        for i in range(len(array) - 1, 0 - 1, -1):
            pattern_input: Tuple[_RESULT, _STATE] = mapping(array[i], acc)
            res[i] = pattern_input[0]
            acc = pattern_input[1]
        return (res, acc)



def indexed(source: Array[_T]) -> Array[Tuple[int, _T]]:
    len_1: int = len(source) or 0
    target: Array[Tuple[int, _T]] = [None]*len_1
    for i in range(0, (len_1 - 1) + 1, 1):
        target[i] = (i, source[i])
    return target


def truncate(count: int, array: Array[_T]) -> Array[_T]:
    def _arrow43(x: int, y: int, count: int=count, array: Array[Any]=array) -> int:
        return compare_primitives(x, y)

    count_1: int = max_1(_arrow43, 0, count) or 0
    return array[0:0+count_1]


def concat(arrays: IEnumerable_1[Array[_T]], cons: Cons_1[_T]) -> Array[_T]:
    arrays_1: Array[Array[_T]] = arrays if (isinstance(arrays, list)) else (list(arrays))
    match_value: int = len(arrays_1) or 0
    if match_value == 0:
        return Helpers_allocateArrayFromCons(cons, 0)

    elif match_value == 1:
        return arrays_1[0]

    else: 
        total_idx: int = 0
        total_length: int = 0
        for idx in range(0, (len(arrays_1) - 1) + 1, 1):
            arr_1: Array[_T] = arrays_1[idx]
            total_length = (total_length + len(arr_1)) or 0
        result: Array[_T] = Helpers_allocateArrayFromCons(cons, total_length)
        for idx_1 in range(0, (len(arrays_1) - 1) + 1, 1):
            arr_2: Array[_T] = arrays_1[idx_1]
            for j in range(0, (len(arr_2) - 1) + 1, 1):
                result[total_idx] = arr_2[j]
                total_idx = (total_idx + 1) or 0
        return result



def collect(mapping: Callable[[_T], Array[_U]], array: Array[_T], cons: Cons_1[_U]) -> Array[_U]:
    return concat(map(mapping, array, None), cons)


def where(predicate: Callable[[__A], bool], array: Array[__A]) -> Array[__A]:
    return list(builtins.filter(predicate, array))


def index_of(array: Array[_T], item_1: _T, start: Optional[int], count: Optional[int], eq: IEqualityComparer_1[Any]) -> int:
    start_1: int = default_arg(start, 0) or 0
    def mapping(c: int, array: Array[Any]=array, item_1: _T=item_1, start: Any=start, count: Any=count, eq: Any=eq) -> int:
        return start_1 + c

    end_0027: int = default_arg(map_1(mapping, count), len(array)) or 0
    def loop(i_mut: int, array: Array[Any]=array, item_1: _T=item_1, start: Any=start, count: Any=count, eq: Any=eq) -> int:
        while True:
            (i,) = (i_mut,)
            if i >= end_0027:
                return -1

            elif eq.Equals(item_1, array[i]):
                return i

            else: 
                i_mut = i + 1
                continue

            break

    return loop(start_1)


def contains(value: _T, array: Array[_T], eq: IEqualityComparer_1[Any]) -> bool:
    return index_of(array, value, None, None, eq) >= 0


def empty(cons: Cons_1[__A]) -> Array[__A]:
    return Helpers_allocateArrayFromCons(cons, 0)


def singleton(value: _T, cons: Cons_1[_T]) -> Array[_T]:
    ar: Array[_T] = Helpers_allocateArrayFromCons(cons, 1)
    ar[0] = value
    return ar


def initialize(count: int, initializer: Callable[[int], _T], cons: Cons_1[_T]) -> Array[_T]:
    if count < 0:
        raise Exception("The input must be non-negative\\nParameter name: count")

    result: Array[_T] = Helpers_allocateArrayFromCons(cons, count)
    for i in range(0, (count - 1) + 1, 1):
        result[i] = initializer(i)
    return result


def pairwise(array: Array[_T]) -> Array[Tuple[_T, _T]]:
    if len(array) < 2:
        return []

    else: 
        count: int = (len(array) - 1) or 0
        result: Array[Tuple[_T, _T]] = [None]*count
        for i in range(0, (count - 1) + 1, 1):
            result[i] = (array[i], array[i + 1])
        return result



def replicate(count: int, initial: _T, cons: Cons_1[_T]) -> Array[_T]:
    if count < 0:
        raise Exception("The input must be non-negative\\nParameter name: count")

    result: Array[_T] = Helpers_allocateArrayFromCons(cons, count)
    for i in range(0, (len(result) - 1) + 1, 1):
        result[i] = initial
    return result


def copy(array: Array[_T]) -> Array[_T]:
    return array[:]


def copy_to(source: Array[_T], source_index: int, target: Array[_T], target_index: int, count: int) -> None:
    copy_to_array(source, source_index, target, target_index, count)


def reverse(array: Array[_T]) -> Array[_T]:
    array_1: Array[_T] = array[:]
    return array_1[::-1]


def scan(folder: Callable[[_STATE, _T], _STATE], state: _STATE, array: Array[_T], cons: Cons_1[_STATE]) -> Array[_STATE]:
    res: Array[_STATE] = Helpers_allocateArrayFromCons(cons, len(array) + 1)
    res[0] = state
    for i in range(0, (len(array) - 1) + 1, 1):
        res[i + 1] = folder(res[i], array[i])
    return res


def scan_back(folder: Callable[[_T, _STATE], _STATE], array: Array[_T], state: _STATE, cons: Cons_1[_STATE]) -> Array[_STATE]:
    res: Array[_STATE] = Helpers_allocateArrayFromCons(cons, len(array) + 1)
    res[len(array)] = state
    for i in range(len(array) - 1, 0 - 1, -1):
        res[i] = folder(array[i], res[i + 1])
    return res


def skip(count: int, array: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    if count > len(array):
        raise Exception("count is greater than array length\\nParameter name: count")

    if count == len(array):
        return Helpers_allocateArrayFromCons(cons, 0)

    else: 
        count_1: int = (0 if (count < 0) else count) or 0
        return array[count_1:]



def skip_while(predicate: Callable[[_T], bool], array: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    count: int = 0
    while predicate(array[count]) if (count < len(array)) else False:
        count = (count + 1) or 0
    if count == len(array):
        return Helpers_allocateArrayFromCons(cons, 0)

    else: 
        return array[count:]



def take(count: int, array: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    if count < 0:
        raise Exception("The input must be non-negative\\nParameter name: count")

    if count > len(array):
        raise Exception("count is greater than array length\\nParameter name: count")

    if count == 0:
        return Helpers_allocateArrayFromCons(cons, 0)

    else: 
        return array[0:0+count]



def take_while(predicate: Callable[[_T], bool], array: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    count: int = 0
    while predicate(array[count]) if (count < len(array)) else False:
        count = (count + 1) or 0
    if count == 0:
        return Helpers_allocateArrayFromCons(cons, 0)

    else: 
        return array[0:0+count]



def add_in_place(x: _T, array: Array[_T]) -> None:
    ignore(array.append(x))


def add_range_in_place(range: IEnumerable_1[_T], array: Array[_T]) -> None:
    with get_enumerator(range) as enumerator:
        while enumerator.System_Collections_IEnumerator_MoveNext():
            add_in_place(enumerator.System_Collections_Generic_IEnumerator_1_get_Current(), array)


def insert_range_in_place(index: int, range: IEnumerable_1[_T], array: Array[_T]) -> None:
    i: int = index or 0
    with get_enumerator(range) as enumerator:
        while enumerator.System_Collections_IEnumerator_MoveNext():
            x: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
            def _arrow44(__unit: None=None, index: int=index, range: IEnumerable_1[Any]=range, array: Array[Any]=array) -> Array[_T]:
                index_1: int = i or 0
                return array.insert(index_1, x)

            ignore(_arrow44())
            i = (i + 1) or 0


def remove_in_place(item_1: _T, array: Array[_T], eq: IEqualityComparer_1[Any]) -> bool:
    i: int = index_of(array, item_1, None, None, eq) or 0
    if i > -1:
        ignore(Helpers_spliceImpl(array, i, 1))
        return True

    else: 
        return False



def remove_all_in_place(predicate: Callable[[_T], bool], array: Array[_T]) -> int:
    def count_remove_all(count: int, predicate: Any=predicate, array: Array[Any]=array) -> int:
        i: int = (next((i for i, x in enumerate(array) if (predicate)(x)), -1)) or 0
        if i > -1:
            ignore(Helpers_spliceImpl(array, i, 1))
            return count_remove_all(count) + 1

        else: 
            return count


    return count_remove_all(0)


def partition(f: Callable[[_T], bool], source: Array[_T], cons: Cons_1[_T]) -> Tuple[Array[_T], Array[_T]]:
    len_1: int = len(source) or 0
    res1: Array[_T] = Helpers_allocateArrayFromCons(cons, len_1)
    res2: Array[_T] = Helpers_allocateArrayFromCons(cons, len_1)
    i_true: int = 0
    i_false: int = 0
    for i in range(0, (len_1 - 1) + 1, 1):
        if f(source[i]):
            res1[i_true] = source[i]
            i_true = (i_true + 1) or 0

        else: 
            res2[i_false] = source[i]
            i_false = (i_false + 1) or 0

    return (truncate(i_true, res1), truncate(i_false, res2))


def find(predicate: Callable[[_T], bool], array: Array[_T]) -> _T:
    match_value: Optional[_T] = next((x for x in array if (predicate)(x)), None)
    if match_value is None:
        return index_not_found()

    else: 
        return value_1(match_value)



def try_find(predicate: Callable[[_T], bool], array: Array[_T]) -> Optional[_T]:
    return next((x for x in array if (predicate)(x)), None)


def find_index(predicate: Callable[[_T], bool], array: Array[_T]) -> int:
    match_value: int = (next((i for i, x in enumerate(array) if (predicate)(x)), -1)) or 0
    if match_value > -1:
        return match_value

    else: 
        index_not_found()
        return -1



def try_find_index(predicate: Callable[[_T], bool], array: Array[_T]) -> Optional[int]:
    match_value: int = (next((i for i, x in enumerate(array) if (predicate)(x)), -1)) or 0
    if match_value > -1:
        return match_value

    else: 
        return None



def pick(chooser: Callable[[__A], Optional[__B]], array: Array[__A]) -> __B:
    def loop(i_mut: int, chooser: Any=chooser, array: Array[Any]=array) -> __B:
        while True:
            (i,) = (i_mut,)
            if i >= len(array):
                return index_not_found()

            else: 
                match_value: Optional[__B] = chooser(array[i])
                if match_value is not None:
                    return value_1(match_value)

                else: 
                    i_mut = i + 1
                    continue


            break

    return loop(0)


def try_pick(chooser: Callable[[__A], Optional[__B]], array: Array[__A]) -> Optional[__B]:
    def loop(i_mut: int, chooser: Any=chooser, array: Array[Any]=array) -> Optional[__B]:
        while True:
            (i,) = (i_mut,)
            if i >= len(array):
                return None

            else: 
                match_value: Optional[__B] = chooser(array[i])
                if match_value is None:
                    i_mut = i + 1
                    continue

                else: 
                    return match_value


            break

    return loop(0)


def find_back(predicate: Callable[[__A], bool], array: Array[__A]) -> __A:
    def loop(i_mut: int, predicate: Any=predicate, array: Array[Any]=array) -> __A:
        while True:
            (i,) = (i_mut,)
            if i < 0:
                return index_not_found()

            elif predicate(array[i]):
                return array[i]

            else: 
                i_mut = i - 1
                continue

            break

    return loop(len(array) - 1)


def try_find_back(predicate: Callable[[__A], bool], array: Array[__A]) -> Optional[__A]:
    def loop(i_mut: int, predicate: Any=predicate, array: Array[Any]=array) -> Optional[__A]:
        while True:
            (i,) = (i_mut,)
            if i < 0:
                return None

            elif predicate(array[i]):
                return some(array[i])

            else: 
                i_mut = i - 1
                continue

            break

    return loop(len(array) - 1)


def find_last_index(predicate: Callable[[__A], bool], array: Array[__A]) -> int:
    def loop(i_mut: int, predicate: Any=predicate, array: Array[Any]=array) -> int:
        while True:
            (i,) = (i_mut,)
            if i < 0:
                return -1

            elif predicate(array[i]):
                return i

            else: 
                i_mut = i - 1
                continue

            break

    return loop(len(array) - 1)


def find_index_back(predicate: Callable[[__A], bool], array: Array[__A]) -> int:
    def loop(i_mut: int, predicate: Any=predicate, array: Array[Any]=array) -> int:
        while True:
            (i,) = (i_mut,)
            if i < 0:
                index_not_found()
                return -1

            elif predicate(array[i]):
                return i

            else: 
                i_mut = i - 1
                continue

            break

    return loop(len(array) - 1)


def try_find_index_back(predicate: Callable[[__A], bool], array: Array[__A]) -> Optional[int]:
    def loop(i_mut: int, predicate: Any=predicate, array: Array[Any]=array) -> Optional[int]:
        while True:
            (i,) = (i_mut,)
            if i < 0:
                return None

            elif predicate(array[i]):
                return i

            else: 
                i_mut = i - 1
                continue

            break

    return loop(len(array) - 1)


def choose(chooser: Callable[[_T], Optional[_U]], array: Array[_T], cons: Cons_1[_U]) -> Array[_U]:
    res: Array[_U] = []
    for i in range(0, (len(array) - 1) + 1, 1):
        match_value: Optional[_U] = chooser(array[i])
        if match_value is not None:
            y: _U = value_1(match_value)
            ignore(res.append(y))

    if equals_1(cons, None):
        return res

    else: 
        def _arrow45(x: Optional[_U]=None, chooser: Any=chooser, array: Array[Any]=array, cons: Any=cons) -> Optional[_U]:
            return x

        return map(_arrow45, res, cons)



def fold_indexed(folder: Callable[[int, _STATE, _T], _STATE], state: _STATE, array: Array[_T]) -> _STATE:
    def _arrow46(delegate_arg: _STATE, delegate_arg_1: _T, delegate_arg_2: int, folder: Any=folder, state: _STATE=state, array: Array[Any]=array) -> _STATE:
        return folder(delegate_arg_2, delegate_arg, delegate_arg_1)

    return array.reduce((_arrow46, state))


def fold(folder: Callable[[_STATE, _T], _STATE], state: _STATE, array: Array[_T]) -> _STATE:
    def _arrow47(acc: _STATE, x: _T, folder: Any=folder, state: _STATE=state, array: Array[Any]=array) -> _STATE:
        return folder(acc, x)

    return functools.reduce(_arrow47, array, state)


def iterate(action: Callable[[_T], None], array: Array[_T]) -> None:
    for i in range(0, (len(array) - 1) + 1, 1):
        action(array[i])


def iterate_indexed(action: Callable[[int, _T], None], array: Array[_T]) -> None:
    for i in range(0, (len(array) - 1) + 1, 1):
        action(i, array[i])


def iterate2(action: Callable[[_T, _T], None], array1: Array[_T], array2: Array[_T]) -> None:
    if len(array1) != len(array2):
        different_lengths()

    for i in range(0, (len(array1) - 1) + 1, 1):
        action(array1[i], array2[i])


def iterate_indexed2(action: Callable[[int, _T, _T], None], array1: Array[_T], array2: Array[_T]) -> None:
    if len(array1) != len(array2):
        different_lengths()

    for i in range(0, (len(array1) - 1) + 1, 1):
        action(i, array1[i], array2[i])


def is_empty(array: Array[Any]) -> bool:
    return len(array) == 0


def for_all(predicate: Callable[[_T], bool], array: Array[_T]) -> bool:
    return all([predicate(x) for x in array])


def permute(f: Callable[[int], int], array: Array[_T]) -> Array[_T]:
    size: int = len(array) or 0
    res: Array[_T] = array[:]
    check_flags: Array[int] = [None]*size
    def _arrow48(i: int, x: _T, f: Any=f, array: Array[Any]=array) -> None:
        j: int = f(i) or 0
        if True if (j < 0) else (j >= size):
            raise Exception("Not a valid permutation")

        res[j] = x
        check_flags[j] = 1

    iterate_indexed(_arrow48, array)
    def predicate(y: int, f: Any=f, array: Array[Any]=array) -> bool:
        return 1 == y

    if not (all([predicate(x) for x in check_flags])):
        raise Exception("Not a valid permutation")

    return res


def set_slice(target: Array[_T], lower: Optional[int], upper: Optional[int], source: Array[_T]) -> None:
    lower_1: int = default_arg(lower, 0) or 0
    upper_1: int = default_arg(upper, -1) or 0
    length: int = ((upper_1 if (upper_1 >= 0) else (len(target) - 1)) - lower_1) or 0
    for i in range(0, length + 1, 1):
        target[i + lower_1] = source[i]


def sort_in_place_by(projection: Callable[[_A], _B], xs: Array[_A], comparer: IComparer_1[_B]) -> None:
    def _arrow49(x: _A, y: _A, projection: Any=projection, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(projection(x), projection(y))

    xs.sort()


def sort_in_place(xs: Array[_T], comparer: IComparer_1[_T]) -> None:
    def _arrow50(x: _T, y: _T, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(x, y)

    xs.sort()


def sort(xs: Array[_T], comparer: IComparer_1[_T]) -> Array[_T]:
    xs_1: Array[_T] = xs[:]
    def comparer_1(x: _T, y: _T, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(x, y)

    xs_1.sort()
    return xs_1


def sort_by(projection: Callable[[_A], _B], xs: Array[_A], comparer: IComparer_1[_B]) -> Array[_A]:
    xs_1: Array[_A] = xs[:]
    def comparer_1(x: _A, y: _A, projection: Any=projection, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(projection(x), projection(y))

    xs_1.sort()
    return xs_1


def sort_descending(xs: Array[_T], comparer: IComparer_1[_T]) -> Array[_T]:
    xs_1: Array[_T] = xs[:]
    def comparer_1(x: _T, y: _T, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(x, y) * -1

    xs_1.sort()
    return xs_1


def sort_by_descending(projection: Callable[[_A], _B], xs: Array[_A], comparer: IComparer_1[_B]) -> Array[_A]:
    xs_1: Array[_A] = xs[:]
    def comparer_1(x: _A, y: _A, projection: Any=projection, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> int:
        return comparer.Compare(projection(x), projection(y)) * -1

    xs_1.sort()
    return xs_1


def sort_with(comparer: Callable[[_T, _T], int], xs: Array[_T]) -> Array[_T]:
    comparer_1: Callable[[_T, _T], int] = comparer
    xs_1: Array[_T] = xs[:]
    xs_1.sort()
    return xs_1


def all_pairs(xs: Array[_T1], ys: Array[_T2]) -> Array[Tuple[_T1, _T2]]:
    len1: int = len(xs) or 0
    len2: int = len(ys) or 0
    res: Array[Tuple[_T1, _T2]] = [None]*(len1 * len2)
    for i in range(0, (len(xs) - 1) + 1, 1):
        for j in range(0, (len(ys) - 1) + 1, 1):
            res[(i * len2) + j] = (xs[i], ys[j])
    return res


def unfold(generator: Callable[[_STATE], Optional[Tuple[_T, _STATE]]], state: _STATE) -> Array[_T]:
    res: Array[_T] = []
    def loop(state_1_mut: Optional[_STATE]=None, generator: Any=generator, state: _STATE=state) -> None:
        while True:
            (state_1,) = (state_1_mut,)
            match_value: Optional[Tuple[_T, _STATE]] = generator(state_1)
            if match_value is not None:
                x: _T = match_value[0]
                s: _STATE = match_value[1]
                ignore(res.append(x))
                state_1_mut = s
                continue

            break

    loop(state)
    return res


def unzip(array: Array[Tuple[__A, __B]]) -> Tuple[Array[__A], Array[__B]]:
    len_1: int = len(array) or 0
    res1: Array[__A] = [None]*len_1
    res2: Array[__B] = [None]*len_1
    def _arrow51(i: int, tupled_arg: Tuple[__A, __B], array: Any=array) -> None:
        res1[i] = tupled_arg[0]
        res2[i] = tupled_arg[1]

    iterate_indexed(_arrow51, array)
    return (res1, res2)


def unzip3(array: Array[Tuple[__A, __B, __C]]) -> Tuple[Array[__A], Array[__B], Array[__C]]:
    len_1: int = len(array) or 0
    res1: Array[__A] = [None]*len_1
    res2: Array[__B] = [None]*len_1
    res3: Array[__C] = [None]*len_1
    def _arrow52(i: int, tupled_arg: Tuple[__A, __B, __C], array: Any=array) -> None:
        res1[i] = tupled_arg[0]
        res2[i] = tupled_arg[1]
        res3[i] = tupled_arg[2]

    iterate_indexed(_arrow52, array)
    return (res1, res2, res3)


def zip(array1: Array[_T], array2: Array[_U]) -> Array[Tuple[_T, _U]]:
    if len(array1) != len(array2):
        different_lengths()

    result: Array[Tuple[_T, _U]] = [None]*len(array1)
    for i in range(0, (len(array1) - 1) + 1, 1):
        result[i] = (array1[i], array2[i])
    return result


def zip3(array1: Array[_T], array2: Array[_U], array3: Array[_U]) -> Array[Tuple[_T, _U, _U]]:
    if True if (len(array1) != len(array2)) else (len(array2) != len(array3)):
        different_lengths()

    result: Array[Tuple[_T, _U, _U]] = [None]*len(array1)
    for i in range(0, (len(array1) - 1) + 1, 1):
        result[i] = (array1[i], array2[i], array3[i])
    return result


def chunk_by_size(chunk_size: int, array: Array[_T]) -> Array[Array[_T]]:
    if chunk_size < 1:
        raise Exception("The input must be positive.\\nParameter name: size")

    if len(array) == 0:
        return [[]]

    else: 
        result: Array[Array[_T]] = []
        for x in range(0, (int(ceil(divide(len(array), chunk_size))) - 1) + 1, 1):
            start: int = (x * chunk_size) or 0
            slice: Array[_T] = array[start:start+chunk_size]
            ignore(result.append(slice))
        return result



def split_at(index: int, array: Array[_T]) -> Tuple[Array[_T], Array[_T]]:
    if True if (index < 0) else (index > len(array)):
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    return (array[0:0+index], array[index:])


def compare_with(comparer: Callable[[_T, _T], int], source1: Array[_T], source2: Array[_T]) -> int:
    if source1 is None:
        if source2 is None:
            return 0

        else: 
            return -1


    elif source2 is None:
        return 1

    else: 
        len1: int = len(source1) or 0
        len2: int = len(source2) or 0
        len_1: int = (len1 if (len1 < len2) else len2) or 0
        i: int = 0
        res: int = 0
        while (i < len_1) if (res == 0) else False:
            res = comparer(source1[i], source2[i]) or 0
            i = (i + 1) or 0
        if res != 0:
            return res

        elif len1 > len2:
            return 1

        elif len1 < len2:
            return -1

        else: 
            return 0




def compare_to(comparer: Callable[[_T, _T], int], source1: Array[_T], source2: Array[_T]) -> int:
    if source1 is None:
        if source2 is None:
            return 0

        else: 
            return -1


    elif source2 is None:
        return 1

    else: 
        len1: int = len(source1) or 0
        len2: int = len(source2) or 0
        if len1 > len2:
            return 1

        elif len1 < len2:
            return -1

        else: 
            i: int = 0
            res: int = 0
            while (i < len1) if (res == 0) else False:
                res = comparer(source1[i], source2[i]) or 0
                i = (i + 1) or 0
            return res




def equals_with(equals: Callable[[_T, _T], bool], array1: Array[_T], array2: Array[_T]) -> bool:
    if array1 is None:
        if array2 is None:
            return True

        else: 
            return False


    elif array2 is None:
        return False

    else: 
        i: int = 0
        result: bool = True
        length1: int = len(array1) or 0
        length2: int = len(array2) or 0
        if length1 > length2:
            return False

        elif length1 < length2:
            return False

        else: 
            while result if (i < length1) else False:
                result = equals(array1[i], array2[i])
                i = (i + 1) or 0
            return result




def exactly_one(array: Array[_T]) -> _T:
    if len(array) == 1:
        return array[0]

    elif len(array) == 0:
        raise Exception("The input sequence was empty\\nParameter name: array")

    else: 
        raise Exception("Input array too long\\nParameter name: array")



def try_exactly_one(array: Array[_T]) -> Optional[_T]:
    if len(array) == 1:
        return some(array[0])

    else: 
        return None



def head(array: Array[_T]) -> _T:
    if len(array) == 0:
        raise Exception("The input array was empty\\nParameter name: array")

    else: 
        return array[0]



def try_head(array: Array[_T]) -> Optional[_T]:
    if len(array) == 0:
        return None

    else: 
        return some(array[0])



def tail(array: Array[_T]) -> Array[_T]:
    if len(array) == 0:
        raise Exception("Not enough elements\\nParameter name: array")

    return array[1:]


def item(index: int, array: Array[__A]) -> __A:
    return array[index]


def try_item(index: int, array: Array[_T]) -> Optional[_T]:
    if True if (index < 0) else (index >= len(array)):
        return None

    else: 
        return some(array[index])



def fold_back_indexed(folder: Callable[[int, _T, _STATE], _STATE], array: Array[_T], state: _STATE) -> _STATE:
    def _arrow53(delegate_arg: _STATE, delegate_arg_1: _T, delegate_arg_2: int, folder: Any=folder, array: Array[Any]=array, state: _STATE=state) -> _STATE:
        return folder(delegate_arg_2, delegate_arg_1, delegate_arg)

    return array.reduceRight((_arrow53, state))


def fold_back(folder: Callable[[_T, _STATE], _STATE], array: Array[_T], state: _STATE) -> _STATE:
    def _arrow54(acc: _STATE, x: _T, folder: Any=folder, array: Array[Any]=array, state: _STATE=state) -> _STATE:
        return folder(x, acc)

    return functools.reduce(_arrow54, array[::-1], state)


def fold_indexed2(folder: Callable[[int, __A, __B, __C], __A], state: __A, array1: Array[__B], array2: Array[__C]) -> __A:
    acc: __A = state
    if len(array1) != len(array2):
        raise Exception("Arrays have different lengths")

    for i in range(0, (len(array1) - 1) + 1, 1):
        acc = folder(i, acc, array1[i], array2[i])
    return acc


def fold2(folder: Callable[[_STATE, _T1, _T2], _STATE], state: _STATE, array1: Array[_T1], array2: Array[_T2]) -> _STATE:
    def _arrow55(_arg: int, acc: _STATE, x: _T1, y: _T2, folder: Any=folder, state: _STATE=state, array1: Array[Any]=array1, array2: Array[Any]=array2) -> _STATE:
        return folder(acc, x, y)

    return fold_indexed2(_arrow55, state, array1, array2)


def fold_back_indexed2(folder: Callable[[int, _T1, _T2, _STATE], _STATE], array1: Array[_T1], array2: Array[_T2], state: _STATE) -> _STATE:
    acc: _STATE = state
    if len(array1) != len(array2):
        different_lengths()

    size: int = len(array1) or 0
    for i in range(1, size + 1, 1):
        acc = folder(i - 1, array1[size - i], array2[size - i], acc)
    return acc


def fold_back2(f: Callable[[_T1, _T2, _STATE], _STATE], array1: Array[_T1], array2: Array[_T2], state: _STATE) -> _STATE:
    def _arrow56(_arg: int, x: _T1, y: _T2, acc: _STATE, f: Any=f, array1: Array[Any]=array1, array2: Array[Any]=array2, state: _STATE=state) -> _STATE:
        return f(x, y, acc)

    return fold_back_indexed2(_arrow56, array1, array2, state)


def reduce(reduction: Callable[[_T, _T], _T], array: Array[_T]) -> _T:
    if len(array) == 0:
        raise Exception("The input array was empty")

    return functools.reduce(reduction, array)


def reduce_back(reduction: Callable[[_T, _T], _T], array: Array[_T]) -> _T:
    if len(array) == 0:
        raise Exception("The input array was empty")

    return functools.reduce(reduction, array[::-1])


def for_all2(predicate: Callable[[__A, __B], bool], array1: Array[__A], array2: Array[__B]) -> bool:
    def _arrow57(acc: bool, x: __A, y: __B, predicate: Any=predicate, array1: Array[Any]=array1, array2: Array[Any]=array2) -> bool:
        return predicate(x, y) if acc else False

    return fold2(_arrow57, True, array1, array2)


def exists_offset(predicate_mut: Callable[[_T], bool], array_mut: Array[_T], index_mut: int) -> bool:
    while True:
        (predicate, array, index) = (predicate_mut, array_mut, index_mut)
        if index == len(array):
            return False

        elif predicate(array[index]):
            return True

        else: 
            predicate_mut = predicate
            array_mut = array
            index_mut = index + 1
            continue

        break


def exists(predicate: Callable[[__A], bool], array: Array[__A]) -> bool:
    return exists_offset(predicate, array, 0)


def exists_offset2(predicate_mut: Callable[[__A, __B], bool], array1_mut: Array[__A], array2_mut: Array[__B], index_mut: int) -> bool:
    while True:
        (predicate, array1, array2, index) = (predicate_mut, array1_mut, array2_mut, index_mut)
        if index == len(array1):
            return False

        elif predicate(array1[index], array2[index]):
            return True

        else: 
            predicate_mut = predicate
            array1_mut = array1
            array2_mut = array2
            index_mut = index + 1
            continue

        break


def exists2(predicate: Callable[[__A, __B], bool], array1: Array[__A], array2: Array[__B]) -> bool:
    if len(array1) != len(array2):
        different_lengths()

    return exists_offset2(predicate, array1, array2, 0)


def sum(array: Array[_T], adder: IGenericAdder_1[_T]) -> _T:
    acc: _T = adder.GetZero()
    for i in range(0, (len(array) - 1) + 1, 1):
        acc = adder.Add(acc, array[i])
    return acc


def sum_by(projection: Callable[[_T], _T2], array: Array[_T], adder: IGenericAdder_1[_T2]) -> _T2:
    acc: _T2 = adder.GetZero()
    for i in range(0, (len(array) - 1) + 1, 1):
        acc = adder.Add(acc, projection(array[i]))
    return acc


def max_by(projection: Callable[[_A], _B], xs: Array[_A], comparer: IComparer_1[_B]) -> _A:
    def _arrow58(x: _A, y: _A, projection: Any=projection, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> _A:
        return y if (comparer.Compare(projection(y), projection(x)) > 0) else x

    return reduce(_arrow58, xs)


def max(xs: Array[_A], comparer: IComparer_1[_A]) -> _A:
    def _arrow59(x: _A, y: _A, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> _A:
        return y if (comparer.Compare(y, x) > 0) else x

    return reduce(_arrow59, xs)


def min_by(projection: Callable[[_A], _B], xs: Array[_A], comparer: IComparer_1[_B]) -> _A:
    def _arrow60(x: _A, y: _A, projection: Any=projection, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> _A:
        return x if (comparer.Compare(projection(y), projection(x)) > 0) else y

    return reduce(_arrow60, xs)


def min(xs: Array[_A], comparer: IComparer_1[_A]) -> _A:
    def _arrow61(x: _A, y: _A, xs: Array[Any]=xs, comparer: IComparer_1[Any]=comparer) -> _A:
        return x if (comparer.Compare(y, x) > 0) else y

    return reduce(_arrow61, xs)


def average(array: Array[_T], averager: IGenericAverager_1[_T]) -> _T:
    if len(array) == 0:
        raise Exception("The input array was empty\\nParameter name: array")

    total: _T = averager.GetZero()
    for i in range(0, (len(array) - 1) + 1, 1):
        total = averager.Add(total, array[i])
    return averager.DivideByInt(total, len(array))


def average_by(projection: Callable[[_T], _T2], array: Array[_T], averager: IGenericAverager_1[_T2]) -> _T2:
    if len(array) == 0:
        raise Exception("The input array was empty\\nParameter name: array")

    total: _T2 = averager.GetZero()
    for i in range(0, (len(array) - 1) + 1, 1):
        total = averager.Add(total, projection(array[i]))
    return averager.DivideByInt(total, len(array))


def windowed(window_size: int, source: Array[_T]) -> Array[Array[_T]]:
    if window_size <= 0:
        raise Exception("windowSize must be positive")

    res: Array[Array[_T]]
    def _arrow62(x: int, y: int, window_size: int=window_size, source: Array[Any]=source) -> int:
        return compare_primitives(x, y)

    len_1: int = max_1(_arrow62, 0, (len(source) - window_size) + 1) or 0
    res = [None]*len_1
    for i in range(window_size, len(source) + 1, 1):
        res[i - window_size] = source[i - window_size:(i - 1) + 1]
    return res


def split_into(chunks: int, array: Array[_T]) -> Array[Array[_T]]:
    if chunks < 1:
        raise Exception("The input must be positive.\\nParameter name: chunks")

    if len(array) == 0:
        return [[]]

    else: 
        result: Array[Array[_T]] = []
        def _arrow63(x: int, y: int, chunks: int=chunks, array: Array[Any]=array) -> int:
            return compare_primitives(x, y)

        chunks_1: int = min_1(_arrow63, chunks, len(array)) or 0
        min_chunk_size: int = (len(array) // chunks_1) or 0
        chunks_with_extra_item: int = (len(array) % chunks_1) or 0
        for i in range(0, (chunks_1 - 1) + 1, 1):
            chunk_size: int = ((min_chunk_size + 1) if (i < chunks_with_extra_item) else min_chunk_size) or 0
            def _arrow64(x_1: int, y_1: int, chunks: int=chunks, array: Array[Any]=array) -> int:
                return compare_primitives(x_1, y_1)

            start: int = ((i * min_chunk_size) + min_1(_arrow64, chunks_with_extra_item, i)) or 0
            slice: Array[_T] = array[start:start+chunk_size]
            ignore(result.append(slice))
        return result



def transpose(arrays: IEnumerable_1[Array[_T]], cons: Cons_1[_T]) -> Array[Array[_T]]:
    arrays_1: Array[Array[_T]] = arrays if (isinstance(arrays, list)) else (list(arrays))
    len_1: int = len(arrays_1) or 0
    if len_1 == 0:
        return [None]*0

    else: 
        len_inner: int = len(arrays_1[0]) or 0
        def predicate(a: Array[_T], arrays: Any=arrays, cons: Any=cons) -> bool:
            return len(a) == len_inner

        if not for_all(predicate, arrays_1):
            different_lengths()

        result: Array[Array[_T]] = [None]*len_inner
        for i in range(0, (len_inner - 1) + 1, 1):
            result[i] = Helpers_allocateArrayFromCons(cons, len_1)
            for j in range(0, (len_1 - 1) + 1, 1):
                result[i][j] = arrays_1[j][i]
        return result



def insert_at(index: int, y: _T, xs: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    len_1: int = len(xs) or 0
    if True if (index < 0) else (index > len_1):
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    target: Array[_T] = Helpers_allocateArrayFromCons(cons, len_1 + 1)
    for i in range(0, (index - 1) + 1, 1):
        target[i] = xs[i]
    target[index] = y
    for i_1 in range(index, (len_1 - 1) + 1, 1):
        target[i_1 + 1] = xs[i_1]
    return target


def insert_many_at(index: int, ys: IEnumerable_1[_T], xs: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    len_1: int = len(xs) or 0
    if True if (index < 0) else (index > len_1):
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    ys_1: Array[_T] = list(ys)
    len2: int = len(ys_1) or 0
    target: Array[_T] = Helpers_allocateArrayFromCons(cons, len_1 + len2)
    for i in range(0, (index - 1) + 1, 1):
        target[i] = xs[i]
    for i_1 in range(0, (len2 - 1) + 1, 1):
        target[index + i_1] = ys_1[i_1]
    for i_2 in range(index, (len_1 - 1) + 1, 1):
        target[i_2 + len2] = xs[i_2]
    return target


def remove_at(index: int, xs: Array[_T]) -> Array[_T]:
    if True if (index < 0) else (index >= len(xs)):
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    i: int = -1
    def predicate(_arg: Optional[_T]=None, index: int=index, xs: Array[Any]=xs) -> bool:
        nonlocal i
        i = (i + 1) or 0
        return i != index

    return filter(predicate, xs)


def remove_many_at(index: int, count: int, xs: Array[_T]) -> Array[_T]:
    i: int = -1
    status: int = -1
    def predicate(_arg: Optional[_T]=None, index: int=index, count: int=count, xs: Array[Any]=xs) -> bool:
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


    ys: Array[_T] = filter(predicate, xs)
    status_1: int = (1 if (((i + 1) == (index + count)) if (status == 0) else False) else status) or 0
    if status_1 < 1:
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + ("index" if (status_1 < 0) else "count"))

    return ys


def update_at(index: int, y: _T, xs: Array[_T], cons: Cons_1[_T]) -> Array[_T]:
    len_1: int = len(xs) or 0
    if True if (index < 0) else (index >= len_1):
        raise Exception((SR_indexOutOfBounds + "\\nParameter name: ") + "index")

    target: Array[_T] = Helpers_allocateArrayFromCons(cons, len_1)
    for i in range(0, (len_1 - 1) + 1, 1):
        target[i] = y if (i == index) else xs[i]
    return target


__all__ = ["index_not_found", "different_lengths", "append", "filter", "fill", "get_sub_array", "last", "try_last", "map_indexed", "map", "map_indexed2", "map2", "map_indexed3", "map3", "map_fold", "map_fold_back", "indexed", "truncate", "concat", "collect", "where", "index_of", "contains", "empty", "singleton", "initialize", "pairwise", "replicate", "copy", "copy_to", "reverse", "scan", "scan_back", "skip", "skip_while", "take", "take_while", "add_in_place", "add_range_in_place", "insert_range_in_place", "remove_in_place", "remove_all_in_place", "partition", "find", "try_find", "find_index", "try_find_index", "pick", "try_pick", "find_back", "try_find_back", "find_last_index", "find_index_back", "try_find_index_back", "choose", "fold_indexed", "fold", "iterate", "iterate_indexed", "iterate2", "iterate_indexed2", "is_empty", "for_all", "permute", "set_slice", "sort_in_place_by", "sort_in_place", "sort", "sort_by", "sort_descending", "sort_by_descending", "sort_with", "all_pairs", "unfold", "unzip", "unzip3", "zip", "zip3", "chunk_by_size", "split_at", "compare_with", "compare_to", "equals_with", "exactly_one", "try_exactly_one", "head", "try_head", "tail", "item", "try_item", "fold_back_indexed", "fold_back", "fold_indexed2", "fold2", "fold_back_indexed2", "fold_back2", "reduce", "reduce_back", "for_all2", "exists_offset", "exists", "exists_offset2", "exists2", "sum", "sum_by", "max_by", "max", "min_by", "min", "average", "average_by", "windowed", "split_into", "transpose", "insert_at", "insert_many_at", "remove_at", "remove_many_at", "update_at"]

