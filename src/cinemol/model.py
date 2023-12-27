from dataclasses import dataclass
import typing as ty

from cinemol.convex_hull import quick_hull_2d
from cinemol.shapes import Sphere  
from cinemol.style import Color, Fill, FillStyle
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D

import numpy as np

@dataclass
class Node:
    index: int
    center: np.ndarray
    radius: ty.Optional[float] # Is not displayed if None
    fill_color: Color
    fill_style: FillStyle

@dataclass 
class Edge:
    start: int
    end: int
    radius: float
    color: ty.Optional[Color] = None # Color copied from nearest node if None 

def get_node_polygon_vertices(this: Sphere, intersects_with: ty.List[Sphere], resolution: int) -> np.ndarray:
    points = this.generate_points_on_surface(resolution) # [res * res, 3]

    centers = np.array([node.center for node in intersects_with]) # [num_others, 3]
    radii = np.array([[node.radius for node in intersects_with]]) # [num_others, 1]
    
    # Get the points that are visible.
    visible = np.linalg.norm(points[:, np.newaxis, :] - centers[np.newaxis, :, :], axis=-1) - radii # [num_points, num_others]
    visible = np.all(visible > 0, axis=-1) # [num_points]
    visible = points[visible] # [num_visible, 3]
    visible = visible[:, :2] # discard z, [num_visible, 2]

    inds = quick_hull_2d(visible)
    verts = [visible[ind] for ind in inds]
    return verts

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
        min_x, min_y, max_x, max_y = np.inf, np.inf, -np.inf, -np.inf
        for node in self.nodes:
            min_x = min(min_x, node.center[0] - node.radius)
            min_y = min(min_y, node.center[1] - node.radius)
            max_x = max(max_x, node.center[0] + node.radius)
            max_y = max(max_y, node.center[1] + node.radius)
        
        min_x -= margin
        min_y -= margin
        max_x += margin
        max_y += margin

        width = max_x - min_x
        height = max_y - min_y
        return ViewBox(min_x, min_y, width, height)

    def draw(self, res: int, verb: bool = False) -> str:
        view_box = self.calculate_view_box()
        svg = Svg(view_box, version=1.0, encoding="UTF-8")

        # TODO: implement masks/clipping to polygons
        # TODO: implement glossy and cartoon styles (per node/edge)
        # TODO: implement space-filling, ball-and-stick, tube, and wireframe depictions (per node/edge)
        # TODO: add perspective
        # TODO: add rotation function
        # TODO: improve performance like only getting points that are visible, putting numpy where possible, maybe even numba....
        # TODO: only draw when interseting with something that came before... not just giving all items back, but check if works for cylinder as well...

        nodes = sorted(self.nodes, key=lambda node: node.center[2], reverse=False)
        spheres = [Sphere(node.center, node.radius) for node in nodes]

        objects, styles = [], []
        for i, node in enumerate(nodes):
            
            this = spheres[i]
            others = spheres[:i] # Only calculate intersection with previous nodes

            intersects_with = [
                other for other in others 
                if np.linalg.norm(this.center - other.center) < this.radius + other.radius
            ]
            
            if i == 0 or len(intersects_with) == 0:
                circle = Circle2D(f"atom-{node.index}", node.center[:2], node.radius)
                objects.append(circle)

            else:
                points = get_node_polygon_vertices(spheres[i], spheres[:i], res)
                polygon = Polygon2D(f"atom-{node.index}", points)
                objects.append(polygon)

            style = Fill(f"atom-{node.index}", node.center, node.radius, node.fill_color, node.fill_style)
            styles.append(style)

            if verb:
                padding = len(str(len(nodes)))
                print(f"{i}".zfill(padding), end="\r")

        svg_str = svg.to_svg(styles, objects)

        return svg_str