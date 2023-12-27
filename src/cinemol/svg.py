from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing as ty

from cinemol.style import Fill 

import numpy as np

class Shape2D(ABC):
    @abstractmethod
    def to_svg(self) -> str:
        ...

@dataclass 
class Circle2D(Shape2D):
    reference: str
    center: np.ndarray
    radius: float   

    def to_svg(self) -> str:
        cx, cy, r = self.center[0], self.center[1], self.radius
        return f'<circle class="{self.reference}" cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}"/>'

@dataclass
class Line2D(Shape2D):
    reference: str
    start: np.ndarray
    end: np.ndarray

    def to_svg(self) -> str:
        x1, y1, x2, y2 = self.start[0], self.start[1], self.end[0], self.end[1]
        return f'<line class="{self.reference}" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'

@dataclass
class Polygon2D(Shape2D):
    reference: str 
    points: np.ndarray # Ordered list of points

    def to_svg(self) -> str:
        points = " ".join([f"{p[0]:.3f},{p[1]:.3f}" for p in self.points])
        return f'<polygon class="{self.reference}" points="{points}"/>'

@dataclass
class ViewBox:
    min_x: float 
    min_y: float 
    width: float 
    height: float

    def __str__(self) -> str:
        return f"ViewBox(min_x={self.min_x:.3f}, min_y={self.min_y:.3f}, width={self.width:.3f}, height={self.height:.3f})"

    def to_svg(self) -> str:
        return f"viewBox=\"{self.min_x:.3f} {self.min_y:.3f} {self.width:.3f} {self.height:.3f}\""
    
@dataclass 
class Svg:
    view_box: ViewBox  
    version: float = 1.0
    encoding: str = "UTF-8"

    def header(self) -> str:
        return (
            f"<?xml version=\"{self.version}\" encoding=\"{self.encoding}\"?>\n"
            f"<svg xmlns=\"http://www.w3.org/2000/svg\" {self.view_box.to_svg()}>"
        )
    
    def footer(self) -> str:
        return "</svg>"
    
    def to_svg(self, fills: ty.List[Fill], objects: ty.List[Shape2D]) -> str:
        header = self.header()
        footer = self.footer()
        
        styles, definitions = [], []
        for fill in fills:
            style, definition = fill.to_svg()
            
            if style is not None:
                styles.append(style)
            
            if definition is not None:
                definitions.append(definition)

        styles = "\n".join(styles)
        definitions = "\n".join(definitions)
        objects = "\n".join([object.to_svg() for object in objects])

        return f"{header}\n<defs>\n<style>\n{styles}\n</style>\n{definitions}\n</defs>\n{objects}\n{footer}"