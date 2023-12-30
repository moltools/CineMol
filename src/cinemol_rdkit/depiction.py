"""
Defines available depictions.
"""
from enum import Enum, auto

class Style(Enum):
    SpaceFilling = auto()
    BallAndStick = auto()
    Tube = auto()
    Wireframe = auto()