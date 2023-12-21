from __future__ import annotations
from abc import abstractmethod
import builtins
from typing import (Protocol, Any)
from ...fable_library.list import FSharpList
from ...fable_library.util import (IEnumerable_1, IDisposable)

class TextIOBase(Protocol):
    @abstractmethod
    def readline(self, __size: int) -> str:
        ...

    @abstractmethod
    def readlines(self, __hint: int) -> FSharpList[str]:
        ...

    @abstractmethod
    def tell(self) -> int:
        ...

    @abstractmethod
    def write(self, __s: str) -> int:
        ...

    @abstractmethod
    def writelines(self, __lines: IEnumerable_1[str]) -> None:
        ...


class TextIOWrapper(TextIOBase, IDisposable):
    pass

def print(obj: Any | None=None) -> None:
    builtins.print(obj)


__all__ = ["print"]

