from abc import abstractmethod
from typing import (TypeVar, Protocol, Generic, Any)

_T = TypeVar("_T")

class IGenericAdder_1(Protocol, Generic[_T]):
    @abstractmethod
    def Add(self, __arg0: _T, __arg1: _T) -> _T:
        ...

    @abstractmethod
    def GetZero(self) -> _T:
        ...


class IGenericAverager_1(Protocol, Generic[_T]):
    @abstractmethod
    def Add(self, __arg0: _T, __arg1: _T) -> _T:
        ...

    @abstractmethod
    def DivideByInt(self, __arg0: _T, __arg1: int) -> _T:
        ...

    @abstractmethod
    def GetZero(self) -> _T:
        ...


class Symbol_wellknown(Protocol):
    @property
    @abstractmethod
    def Symbol_toStringTag(self) -> str:
        ...


class IJsonSerializable(Protocol):
    @abstractmethod
    def to_json(self) -> Any:
        ...


SR_indexOutOfBounds: str = "The index was outside the range of elements in the collection."

SR_inputWasEmpty: str = "Collection was empty."

SR_inputMustBeNonNegative: str = "The input must be non-negative."

SR_inputSequenceEmpty: str = "The input sequence was empty."

SR_inputSequenceTooLong: str = "The input sequence contains more than one element."

SR_keyNotFoundAlt: str = "An index satisfying the predicate was not found in the collection."

SR_differentLengths: str = "The collections had different lengths."

SR_notEnoughElements: str = "The input sequence has an insufficient number of elements."

__all__ = ["SR_indexOutOfBounds", "SR_inputWasEmpty", "SR_inputMustBeNonNegative", "SR_inputSequenceEmpty", "SR_inputSequenceTooLong", "SR_keyNotFoundAlt", "SR_differentLengths", "SR_notEnoughElements"]

