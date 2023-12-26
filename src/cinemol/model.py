from dataclasses import dataclass
import typing as ty

from cinemol.geometry import Point2D, Point3D, Circle2D    
from cinemol.style import Color, Style 
from cinemol.svg import ViewBox, Svg 

@dataclass
class Node:
    index: int
    position: Point3D
    radius: ty.Optional[float] # Is not displayed if None
    color: Color

@dataclass 
class Edge:
    start: int
    end: int
    radius: float
    color: ty.Optional[Color] = None # Color copied from nearest node if None 

class Scene:
    def __init__(self, nodes: ty.List[Node] = [], edges: ty.List[Edge] = []) -> None:
        self.nodes = nodes
        self.edges = edges

    def __str__(self) -> str:
        return f"Scene(nodes={len(self.nodes)}, edges={len(self.edges)})"

    def add_node(self, node: Node) -> None:
        self.nodes.append(node)
    
    def add_edge(self, edge: Edge) -> None:
        self.edges.append(edge)

    def calculate_view_box(self, margin: float = 5) -> ViewBox:
        min_x = min([node.position.x for node in self.nodes]) - margin
        min_y = min([node.position.y for node in self.nodes]) - margin
        max_x = max([node.position.x for node in self.nodes]) + margin
        max_y = max([node.position.y for node in self.nodes]) + margin
        width = max_x - min_x
        height = max_y - min_y
        return ViewBox(min_x, min_y, width, height)

    def draw(self) -> str:
        view_box = self.calculate_view_box()
        svg = Svg(view_box, version=1.0, encoding="UTF-8")

        styles = [Style(f"atom-{node.index}", node.color) for node in self.nodes]

        # TODO: implement glossy and cartoon styles (per node/edge)
        # TODO: implement space-filling, ball-and-stick, tube, and wireframe depictions (per node/edge)
        # TODO: implement masks/clipping to polygons
        # TODO: add perspective
        # TODO: add rotation function

        objects = []
        for node in self.nodes:
            center = Point2D(node.position.x, node.position.y)
            circle = Circle2D(f"atom-{node.index}", center, node.radius)
            objects.append(circle)

        svg_str = svg.to_svg(styles, objects)

        return svg_str