from __future__ import annotations
from abc import abstractmethod
from typing import Protocol
from ...fable_library.list import FSharpList
from ...fable_library.types import (Array, int64)

class ITimeSeries(Protocol):
    @abstractmethod
    def plot(self, start: str, end: str, aggregates: Array[str], granularity: str) -> None:
        ...


class ITimeSeriesApi(Protocol):
    @abstractmethod
    def list(self) -> FSharpList[ITimeSeries]:
        ...

    @abstractmethod
    def retrieve(self, id: int64) -> ITimeSeries:
        ...


