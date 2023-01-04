from __future__ import annotations
from typing import (TypeVar, Any, Generic, Optional)
from .array import (fill, copy_to, initialize)
from .double import divide
from .reflection import (TypeInfo, class_type)
from .seq import (delay, enumerate_while, append, singleton, empty, to_array)
from .types import (Array, FSharpRef)
from .util import (compare, IComparer_1, equals, structural_hash, IEqualityComparer_1, get_enumerator, IEnumerable_1, IEnumerator, to_iterator, max, compare_primitives)

_T = TypeVar("_T")

_T_ = TypeVar("_T_")

def _expr4(gen0: TypeInfo) -> TypeInfo:
    return class_type("System.Collections.Generic.Comparer`1", [gen0], Comparer_1)


class Comparer_1(Generic[_T]):
    def __init__(self, __unit: None=None) -> None:
        pass

    def Compare(self, x: _T, y: _T) -> int:
        return compare(x, y)


Comparer_1_reflection = _expr4

def Comparer_1__ctor(__unit: None=None) -> Comparer_1[Any]:
    return Comparer_1(__unit)


def Comparer_1_get_Default(__unit: None=None) -> IComparer_1[Any]:
    class ObjectExpr5(IComparer_1[_T_]):
        def Compare(self, x: _T_, y: _T_) -> int:
            return compare(x, y)

    return ObjectExpr5()


def _expr6(gen0: TypeInfo) -> TypeInfo:
    return class_type("System.Collections.Generic.EqualityComparer`1", [gen0], EqualityComparer_1)


class EqualityComparer_1(Generic[_T]):
    def __init__(self, __unit: None=None) -> None:
        pass

    def __eq__(self, x: _T, y: _T) -> bool:
        return equals(x, y)

    def GetHashCode(self, x: Optional[_T]=None) -> int:
        return structural_hash(x)


EqualityComparer_1_reflection = _expr6

def EqualityComparer_1__ctor(__unit: None=None) -> EqualityComparer_1[Any]:
    return EqualityComparer_1(__unit)


def EqualityComparer_1_get_Default(__unit: None=None) -> IEqualityComparer_1[Any]:
    class ObjectExpr7(IEqualityComparer_1[Any]):
        def Equals(self, x: _T_, y: _T_) -> bool:
            return equals(x, y)

        def GetHashCode(self, x_1: Optional[_T_]=None) -> int:
            return structural_hash(x_1)

    return ObjectExpr7()


def _expr12(gen0: TypeInfo) -> TypeInfo:
    return class_type("System.Collections.Generic.Stack`1", [gen0], Stack_1)


class Stack_1(Generic[_T]):
    def __init__(self, initial_contents: Array[Any], initial_count: int) -> None:
        self.contents: Array[_T] = initial_contents
        self.count: int = initial_count or 0

    def GetEnumerator(self, __unit: None=None) -> IEnumerator[_T]:
        _: Stack_1[_T] = self
        def _arrow11(__unit: None=None) -> IEnumerable_1[_T]:
            index: int = (_.count - 1) or 0
            def _arrow8(__unit: None=None) -> bool:
                return index >= 0

            def _arrow10(__unit: None=None) -> IEnumerable_1[_T]:
                def _arrow9(__unit: None=None) -> IEnumerable_1[_T]:
                    nonlocal index
                    index = (index - 1) or 0
                    return empty()

                return append(singleton(_.contents[index]), delay(_arrow9))

            return enumerate_while(_arrow8, delay(_arrow10))

        return get_enumerator(delay(_arrow11))

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        this: Stack_1[_T] = self
        return get_enumerator(this)


Stack_1_reflection = _expr12

def Stack_1__ctor_Z3B4C077E(initial_contents: Array[Any], initial_count: int) -> Stack_1[_T]:
    return Stack_1(initial_contents, initial_count)


def Stack_1__ctor_Z524259A4(initial_capacity: int) -> Stack_1[Any]:
    return Stack_1__ctor_Z3B4C077E(fill([0] * initial_capacity, 0, initial_capacity, None), 0)


def Stack_1__ctor(__unit: None=None) -> Stack_1[Any]:
    return Stack_1__ctor_Z524259A4(4)


def Stack_1__ctor_BB573A(xs: IEnumerable_1[_T]) -> Stack_1[_T]:
    arr: Array[_T] = list(xs)
    return Stack_1__ctor_Z3B4C077E(arr, len(arr))


def Stack_1__Ensure_Z524259A4(_: Stack_1[Any], new_size: int) -> None:
    old_size: int = len(_.contents) or 0
    if new_size > old_size:
        old: Array[_T] = _.contents
        def _arrow15(x: int, y: int, _: Any=_, new_size: int=new_size) -> int:
            return compare_primitives(x, y)

        def _arrow16(x: int, y: int, _: Any=_, new_size: int=new_size) -> int:
            return compare_primitives(x, y)

        _.contents = fill([0] * max(_arrow15, new_size, old_size * 2), 0, max(_arrow16, new_size, old_size * 2), None)
        copy_to(old, 0, _.contents, 0, _.count)



def Stack_1__get_Count(_: Stack_1[Any]) -> int:
    return _.count


def Stack_1__Pop(_: Stack_1[_T]) -> _T:
    _.count = (_.count - 1) or 0
    return _.contents[_.count]


def Stack_1__Peek(_: Stack_1[_T]) -> _T:
    return _.contents[_.count - 1]


def Stack_1__Contains_2B595(_: Stack_1[_T], x: _T) -> bool:
    found: bool = False
    i: int = 0
    while (not found) if (i < _.count) else False:
        if equals(x, _.contents[i]):
            found = True

        else: 
            i = (i + 1) or 0

    return found


def Stack_1__TryPeek_1F3DB691(this: Stack_1[_T], result: FSharpRef[_T]) -> bool:
    if this.count > 0:
        result.contents = Stack_1__Peek(this)
        return True

    else: 
        return False



def Stack_1__TryPop_1F3DB691(this: Stack_1[_T], result: FSharpRef[_T]) -> bool:
    if this.count > 0:
        result.contents = Stack_1__Pop(this)
        return True

    else: 
        return False



def Stack_1__Push_2B595(this: Stack_1[_T], x: _T) -> None:
    Stack_1__Ensure_Z524259A4(this, this.count + 1)
    this.contents[this.count] = x
    this.count = (this.count + 1) or 0


def Stack_1__Clear(_: Stack_1[Any]) -> None:
    _.count = 0
    fill(_.contents, 0, len(_.contents), None)


def Stack_1__TrimExcess(this: Stack_1[Any]) -> None:
    if divide(this.count, len(this.contents)) > 0.9:
        Stack_1__Ensure_Z524259A4(this, this.count)



def Stack_1__ToArray(_: Stack_1[_T]) -> Array[_T]:
    def _arrow17(i: int, _: Any=_) -> _T:
        return _.contents[(_.count - 1) - i]

    return initialize(_.count, _arrow17, None)


def _expr18(gen0: TypeInfo) -> TypeInfo:
    return class_type("System.Collections.Generic.Queue`1", [gen0], Queue_1)


class Queue_1(Generic[_T]):
    def __init__(self, initial_contents: Array[Any], initial_count: int) -> None:
        self.contents: Array[_T] = initial_contents
        self.count: int = initial_count or 0
        self.head: int = 0
        self.tail: int = initial_count or 0

    def GetEnumerator(self, __unit: None=None) -> IEnumerator[_T]:
        _: Queue_1[_T] = self
        return get_enumerator(Queue_1__toSeq(_))

    def __iter__(self) -> IEnumerator[_T]:
        return to_iterator(self.GetEnumerator())

    def System_Collections_IEnumerable_GetEnumerator(self, __unit: None=None) -> IEnumerator[Any]:
        this: Queue_1[_T] = self
        return get_enumerator(this)


Queue_1_reflection = _expr18

def Queue_1__ctor_Z3B4C077E(initial_contents: Array[Any], initial_count: int) -> Queue_1[_T]:
    return Queue_1(initial_contents, initial_count)


def Queue_1__ctor_Z524259A4(initial_capacity: int) -> Queue_1[Any]:
    if initial_capacity < 0:
        raise Exception("capacity is less than 0")

    return Queue_1__ctor_Z3B4C077E(fill([0] * initial_capacity, 0, initial_capacity, None), 0)


def Queue_1__ctor(__unit: None=None) -> Queue_1[Any]:
    return Queue_1__ctor_Z524259A4(4)


def Queue_1__ctor_BB573A(xs: IEnumerable_1[_T]) -> Queue_1[_T]:
    arr: Array[_T] = list(xs)
    return Queue_1__ctor_Z3B4C077E(arr, len(arr))


def Queue_1__get_Count(_: Queue_1[Any]) -> int:
    return _.count


def Queue_1__Enqueue_2B595(_: Queue_1[_T], value: _T) -> None:
    if _.count == Queue_1__size(_):
        Queue_1__ensure_Z524259A4(_, _.count + 1)

    _.contents[_.tail] = value
    _.tail = ((_.tail + 1) % Queue_1__size(_)) or 0
    _.count = (_.count + 1) or 0


def Queue_1__Dequeue(_: Queue_1[_T]) -> _T:
    if _.count == 0:
        raise Exception("Queue is empty")

    value: _T = _.contents[_.head]
    _.head = ((_.head + 1) % Queue_1__size(_)) or 0
    _.count = (_.count - 1) or 0
    return value


def Queue_1__Peek(_: Queue_1[_T]) -> _T:
    if _.count == 0:
        raise Exception("Queue is empty")

    return _.contents[_.head]


def Queue_1__TryDequeue_1F3DB691(this: Queue_1[_T], result: FSharpRef[_T]) -> bool:
    if this.count == 0:
        return False

    else: 
        result.contents = Queue_1__Dequeue(this)
        return True



def Queue_1__TryPeek_1F3DB691(this: Queue_1[_T], result: FSharpRef[_T]) -> bool:
    if this.count == 0:
        return False

    else: 
        result.contents = Queue_1__Peek(this)
        return True



def Queue_1__Contains_2B595(_: Queue_1[_T], x: _T) -> bool:
    found: bool = False
    i: int = 0
    while (not found) if (i < _.count) else False:
        if equals(x, _.contents[Queue_1__toIndex_Z524259A4(_, i)]):
            found = True

        else: 
            i = (i + 1) or 0

    return found


def Queue_1__Clear(_: Queue_1[Any]) -> None:
    _.count = 0
    _.head = 0
    _.tail = 0
    fill(_.contents, 0, Queue_1__size(_), None)


def Queue_1__TrimExcess(_: Queue_1[Any]) -> None:
    if divide(_.count, len(_.contents)) > 0.9:
        Queue_1__ensure_Z524259A4(_, _.count)



def Queue_1__ToArray(_: Queue_1[_T]) -> Array[_T]:
    return to_array(Queue_1__toSeq(_))


def Queue_1__CopyTo_Z3B4C077E(_: Queue_1[_T], target: Array[_T], start: int) -> None:
    i: int = start or 0
    with get_enumerator(Queue_1__toSeq(_)) as enumerator:
        while enumerator.System_Collections_IEnumerator_MoveNext():
            item: _T = enumerator.System_Collections_Generic_IEnumerator_1_get_Current()
            target[i] = item
            i = (i + 1) or 0


def Queue_1__size(this: Queue_1[Any]) -> int:
    return len(this.contents)


def Queue_1__toIndex_Z524259A4(this: Queue_1[Any], i: int) -> int:
    return (this.head + i) % Queue_1__size(this)


def Queue_1__ensure_Z524259A4(this: Queue_1[Any], required_size: int) -> None:
    new_buffer: Array[_T] = fill([0] * required_size, 0, required_size, None)
    if this.head < this.tail:
        copy_to(this.contents, this.head, new_buffer, 0, this.count)

    else: 
        copy_to(this.contents, this.head, new_buffer, 0, Queue_1__size(this) - this.head)
        copy_to(this.contents, 0, new_buffer, Queue_1__size(this) - this.head, this.tail)

    this.head = 0
    this.tail = this.count or 0
    this.contents = new_buffer


def Queue_1__toSeq(this: Queue_1[_T]) -> IEnumerable_1[_T]:
    def _arrow23(__unit: None=None, this: Any=this) -> IEnumerable_1[_T]:
        i: int = 0
        def _arrow20(__unit: None=None) -> bool:
            return i < this.count

        def _arrow22(__unit: None=None) -> IEnumerable_1[_T]:
            def _arrow21(__unit: None=None) -> IEnumerable_1[_T]:
                nonlocal i
                i = (i + 1) or 0
                return empty()

            return append(singleton(this.contents[Queue_1__toIndex_Z524259A4(this, i)]), delay(_arrow21))

        return enumerate_while(_arrow20, delay(_arrow22))

    return delay(_arrow23)


__all__ = ["Comparer_1_reflection", "Comparer_1_get_Default", "EqualityComparer_1_reflection", "EqualityComparer_1_get_Default", "Stack_1_reflection", "Stack_1__ctor_Z524259A4", "Stack_1__ctor", "Stack_1__ctor_BB573A", "Stack_1__Ensure_Z524259A4", "Stack_1__get_Count", "Stack_1__Pop", "Stack_1__Peek", "Stack_1__Contains_2B595", "Stack_1__TryPeek_1F3DB691", "Stack_1__TryPop_1F3DB691", "Stack_1__Push_2B595", "Stack_1__Clear", "Stack_1__TrimExcess", "Stack_1__ToArray", "Queue_1_reflection", "Queue_1__ctor_Z524259A4", "Queue_1__ctor", "Queue_1__ctor_BB573A", "Queue_1__get_Count", "Queue_1__Enqueue_2B595", "Queue_1__Dequeue", "Queue_1__Peek", "Queue_1__TryDequeue_1F3DB691", "Queue_1__TryPeek_1F3DB691", "Queue_1__Contains_2B595", "Queue_1__Clear", "Queue_1__TrimExcess", "Queue_1__ToArray", "Queue_1__CopyTo_Z3B4C077E", "Queue_1__size", "Queue_1__toIndex_Z524259A4", "Queue_1__ensure_Z524259A4", "Queue_1__toSeq"]

