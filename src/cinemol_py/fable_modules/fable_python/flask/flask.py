from __future__ import annotations
from abc import abstractmethod
from typing import Protocol

class RequestBase(Protocol):
    @property
    @abstractmethod
    def url(self) -> str:
        ...


class Request(RequestBase):
    pass

class Flask(Protocol):
    pass

class FlaskStatic(Protocol):
    @abstractmethod
    def Create(self, name: str, static_url_path: str) -> Flask:
        ...


class IExports(Protocol):
    @property
    @abstractmethod
    def request(self) -> Request:
        ...

    @abstractmethod
    def url_for(self, route: str, filename: str) -> str:
        ...


