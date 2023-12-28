"""
This module contains classes for generating SVG documents.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing as ty

from cinemol.style import Fill 

import numpy as np

class Shape2D(ABC):
    """
    Abstract base class for 2D shapes.
    """
    @abstractmethod
    def to_svg(self) -> str:
        """
        Return the SVG representation of the shape.

        :return: The SVG representation of the shape.
        :rtype: str
        """
        ...

@dataclass 
class Circle2D(Shape2D):
    """
    A 2D circle.

    :param str reference: The reference name of the circle.
    :param np.ndarray center: The center of the circle with shape (2,).
    :param float radius: The radius of the circle.
    """
    reference: str
    center: np.ndarray
    radius: float   

    def to_svg(self) -> str:
        """
        Return the SVG representation of the circle.

        :return: The SVG representation of the circle.
        :rtype: str
        """
        cx, cy, r = self.center[0], self.center[1], self.radius
        return f'<circle class="{self.reference}" cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}"/>'

@dataclass
class Line2D(Shape2D):
    """
    A 2D line.
    
    :param str reference: The reference name of the line.
    :param np.ndarray start: The start point of the line with shape (2,).
    :param np.ndarray end: The end point of the line with shape (2,).
    """
    reference: str
    start: np.ndarray
    end: np.ndarray

    def to_svg(self) -> str:
        """
        Return the SVG representation of the line.
        
        :return: The SVG representation of the line.
        :rtype: str
        """
        x1, y1, x2, y2 = self.start[0], self.start[1], self.end[0], self.end[1]
        return f'<line class="{self.reference}" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'

@dataclass
class Polygon2D(Shape2D):
    """
    A 2D polygon.
    
    :param str reference: The reference name of the polygon.
    :param np.ndarray points: The points of the polygon with shape (N, 2).
    """
    reference: str 
    points: np.ndarray # Ordered list of points

    def to_svg(self) -> str:
        """
        Return the SVG representation of the polygon.
        
        :return: The SVG representation of the polygon.
        :rtype: str
        """
        points = " ".join([f"{p[0]:.3f},{p[1]:.3f}" for p in self.points])
        return f'<polygon class="{self.reference}" points="{points}"/>'

@dataclass
class ViewBox:
    """
    A view box.
    
    :param float min_x: The minimum x coordinate.
    :param float min_y: The minimum y coordinate.
    :param float width: The width of the view box.
    :param float height: The height of the view box.
    """
    min_x: float 
    min_y: float 
    width: float 
    height: float

    def __str__(self) -> str:
        """
        Return the string representation of the view box.
        
        :return: The string representation of the view box.
        :rtype: str
        """
        return f"ViewBox(min_x={self.min_x:.3f}, min_y={self.min_y:.3f}, width={self.width:.3f}, height={self.height:.3f})"

    def to_svg(self) -> str:
        """
        Return the SVG representation of the view box.
        
        :return: The SVG representation of the view box.
        :rtype: str
        """
        return f"viewBox=\"{self.min_x:.3f} {self.min_y:.3f} {self.width:.3f} {self.height:.3f}\""
    
@dataclass 
class Svg:
    """
    An SVG document.
    
    :param ViewBox view_box: The view box of the SVG document.
    :param float version: The version of the SVG document.
    :param str encoding: The encoding of the SVG document.
    """
    view_box: ViewBox  
    version: float = 1.0
    encoding: str = "UTF-8"

    def header(self) -> str:
        """
        Return the header of the SVG document.
        
        :return: The header of the SVG document.
        :rtype: str
        """
        return f"<?xml version=\"{self.version}\" encoding=\"{self.encoding}\"?>\n" +\
               f"<svg xmlns=\"http://www.w3.org/2000/svg\" {self.view_box.to_svg()}>"
    
    def footer(self) -> str:
        """
        Return the footer of the SVG document.
        
        :return: The footer of the SVG document.
        :rtype: str
        """
        return "</svg>"
    
    def to_svg(self, fills: ty.List[Fill], objects: ty.List[Shape2D]) -> str:
        """
        Return the SVG representation of the SVG document.
        
        :param ty.List[Fill] fills: The fills of the SVG document.
        :param ty.List[Shape2D] objects: The objects of the SVG document.
        :return: The SVG representation of the SVG document.
        :rtype: str
        """
        header = self.header()
        footer = self.footer()
        
        styles, definitions = [], []
        for fill in fills:
            style, definition = fill.to_svg()
            
            if style is not None: # Should never be None.
                styles.append(style)
            
                if definition is not None:
                    definitions.append(definition)

        styles = "\n".join(styles)
        definitions = "\n".join(definitions)
        objects = "\n".join([object.to_svg() for object in objects])

        return f"{header}\n<defs>\n<style>\n{styles}\n</style>\n{definitions}\n</defs>\n{objects}\n{footer}"