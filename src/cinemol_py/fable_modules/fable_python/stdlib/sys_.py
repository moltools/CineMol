from abc import abstractmethod
from typing import Protocol

class VersionInfo(Protocol):
    @property
    @abstractmethod
    def major(self) -> int:
        ...

    @property
    @abstractmethod
    def micro(self) -> int:
        ...

    @property
    @abstractmethod
    def minor(self) -> int:
        ...

    @property
    @abstractmethod
    def releaselevel(self) -> str:
        ...

    @property
    @abstractmethod
    def serial(self) -> int:
        ...


