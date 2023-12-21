from __future__ import annotations
from abc import abstractmethod
from typing import (Generic, Protocol, TypeVar)
from .events import AbstractEventLoop
from .futures import (Future_1, Awaitable_1)

_T = TypeVar("_T")

class Task_1(Future_1, Awaitable_1, Generic[_T]):
    pass

class Task(Future_1, Awaitable_1):
    pass

class Coroutine_1(Awaitable_1, Generic[_T]):
    @abstractmethod
    def close(self) -> None:
        ...

    @abstractmethod
    def send(self, __arg0: _T) -> None:
        ...

    @abstractmethod
    def throw(self, __arg0: Exception) -> None:
        ...


class IExports(Protocol):
    @abstractmethod
    def get_event_loop(self) -> AbstractEventLoop:
        ...

    @abstractmethod
    def get_running_loop(self) -> AbstractEventLoop:
        ...


