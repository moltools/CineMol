from __future__ import annotations
from abc import abstractmethod
from typing import Protocol

class IWidget(Protocol):
    @abstractmethod
    def close(self) -> None:
        ...


class IIntSlider(IWidget):
    @property
    @abstractmethod
    def max(self) -> int:
        ...

    @property
    @abstractmethod
    def min(self) -> int:
        ...

    @property
    @abstractmethod
    def value(self) -> int:
        ...


class IFloatSlider(IWidget):
    @property
    @abstractmethod
    def max(self) -> float:
        ...

    @property
    @abstractmethod
    def min(self) -> float:
        ...

    @property
    @abstractmethod
    def value(self) -> float:
        ...


class IExports(Protocol):
    pass

