from __future__ import annotations
from abc import abstractmethod
from typing import (Protocol, TypeVar)
from .futures import Future_1

_T = TypeVar("_T")

class AbstractEventLoop(Protocol):
    @abstractmethod
    def create_future(self) -> Future_1[_T]:
        ...

    @abstractmethod
    def run_forever(self) -> None:
        ...


