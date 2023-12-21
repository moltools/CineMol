from abc import abstractmethod
from typing import (Any, Protocol)

class IDisplay(Protocol):
    @abstractmethod
    def Markdown(self, data: str) -> None:
        ...

    @abstractmethod
    def display(self, value: Any) -> None:
        ...


