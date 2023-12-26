from dataclasses import dataclass
import typing as ty

from cinemol.style import Style 
from cinemol.geometry import Shape2D

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
    
    def to_svg(self, styles: ty.List[Style], objects: ty.List[Shape2D]) -> str:
        header = self.header()
        footer = self.footer()
        
        styles = "\n".join([style.to_svg() for style in styles])
        objects = "\n".join([object.to_svg() for object in objects])

        return f"{header}\n<defs>\n<style>\n{styles}\n</style>\n</defs>\n{objects}\n{footer}"