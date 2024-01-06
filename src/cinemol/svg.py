"""
This module contains classes for generating SVG documents.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import typing as ty

from cinemol.geometry import Point2D
from cinemol.style import Fill, Color

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
    :param Point2D center: The center of the circle with shape (2,).
    :param float radius: The radius of the circle.
    """
    reference: str
    center: Point2D
    radius: float   

    def to_svg(self) -> str:
        """
        Return the SVG representation of the circle.

        :return: The SVG representation of the circle.
        :rtype: str
        """
        cx, cy, r = self.center.x, self.center.y, self.radius
        return f'<circle class="{self.reference}" cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}"/>'

@dataclass
class Line2D(Shape2D):
    """
    A 2D line.
    
    :param str reference: The reference name of the line.
    :param Point2D start: The start point of the line with shape (2,).
    :param Point2D end: The end point of the line with shape (2,).
    """
    reference: str
    start: Point2D
    end: Point2D

    def to_svg(self) -> str:
        """
        Return the SVG representation of the line.
        
        :return: The SVG representation of the line.
        :rtype: str
        """
        x1, y1, x2, y2 = self.start.x, self.start.y, self.end.x, self.end.y
        return f'<line class="{self.reference}" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'

@dataclass
class Polygon2D(Shape2D):
    """
    A 2D polygon.
    
    :param str reference: The reference name of the polygon.
    :param ty.List[Point2D] points: The vertices of the polygon.
    """
    reference: str 
    points: ty.List[Point2D] # Ordered list of points.

    def to_svg(self) -> str:
        """
        Return the SVG representation of the polygon.
        
        :return: The SVG representation of the polygon.
        :rtype: str
        """
        points = " ".join([f"{p.x:.3f},{p.y:.3f}" for p in self.points])
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
    
class Svg:
    """
    An SVG document.
    """
    def __init__(
        self,
        view_box: ViewBox,
        window: ty.Optional[ty.Tuple[float, float]] = None,
        background_color: ty.Optional[Color] = None,
        version: float = 1.0,
        encoding: str = "UTF-8",
        fills: ty.List[Fill] = [],
        objects: ty.List[Shape2D] = []
    ):
        """
        :param ViewBox view_box: The view box of the SVG document.
        :param ty.Optional[ty.Tuple[float, float]] window: The window of the SVG document.
        :param ty.Optional[Color] background_color: The background color of the SVG document.
        :param float version: The version of the SVG document.
        :param str encoding: The encoding of the SVG document.
        :param ty.List[Fill] fills: The fills of the SVG document.
        :param ty.List[Shape2D] objects: The objects of the SVG document.
        """
        self.view_box = view_box
        self.window = window
        self.background_color = background_color
        self.version = version
        self.encoding = encoding
        self.fills = fills
        self.objects = objects

    def header(self) -> str:
        """
        Return the header of the SVG document.
        
        :return: The header of the SVG document.
        :rtype: str
        """
        x = self.view_box.min_x
        y = self.view_box.min_y

        if self.background_color is not None:
            background = f"<rect x=\"{x:.3f}\" y=\"{y:.3f}\" width=\"100%\" height=\"100%\" fill=\"{self.background_color.to_hex()}\"/>"
        else:
            background = ""

        if not self.window:
            return f"<?xml version=\"{self.version}\" encoding=\"{self.encoding}\"?>\n" +\
                f"<svg xmlns=\"http://www.w3.org/2000/svg\" {self.view_box.to_svg()}>" + background
        
        else:
            width, height = self.window
            return f"<?xml version=\"{self.version}\" encoding=\"{self.encoding}\"?>\n" +\
                f"<svg xmlns=\"http://www.w3.org/2000/svg\" {self.view_box.to_svg()} width=\"{width}\" height=\"{height}\">" + background
        
    def footer(self) -> str:
        """
        Return the footer of the SVG document.
        
        :return: The footer of the SVG document.
        :rtype: str
        """
        return "</svg>"
    
    def to_svg(self) -> str:
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
        for fill in self.fills:
            style, definition = fill.to_svg()
            
            if style is not None: # Should never be None.
                styles.append(style)
            
                if definition is not None:
                    definitions.append(definition)

        styles = "\n".join(styles)
        definitions = "\n".join(definitions)
        objects = "\n".join([object.to_svg() for object in self.objects])

        return f"{header}\n<defs>\n<style>\n{styles}\n</style>\n{definitions}\n</defs>\n{objects}\n{footer}"