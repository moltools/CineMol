from __future__ import annotations
from abc import abstractmethod
from collections.abc import Callable
from typing import (Protocol, Generic, TypeVar)

_T = TypeVar("_T")

class Awaitable_1(Protocol, Generic[_T]):
    pass

class Future_1(Awaitable_1, Generic[_T]):
    @abstractmethod
    def add_done_callback(self, callback: Callable[[], None]) -> None:
        ...

    @abstractmethod
    def cancel(self) -> None:
        ...

    @abstractmethod
    def cancelled(self) -> bool:
        ...

    @abstractmethod
    def done(self) -> bool:
        ...

    @abstractmethod
    def exception(self) -> Exception:
        ...

    @abstractmethod
    def get_name(self) -> str:
        ...

    @abstractmethod
    def result(self) -> _T:
        ...

    @abstractmethod
    def set_exception(self, __exception: Exception) -> None:
        ...

    @abstractmethod
    def set_name(self, __arg0: str) -> None:
        ...

    @abstractmethod
    def set_result(self, __result: _T) -> None:
        ...


