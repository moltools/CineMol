"""
Contains definitions for a model scene.
"""
from abc import ABC, abstractmethod
import typing as ty

from cinemol.geometry import Point2D
from cinemol.fitting import calculate_convex_hull
from cinemol.shapes import Line3D, Sphere, Cylinder 
from cinemol.style import Color, Fill, FillStyleType, Solid, RadialGradient, LinearGradient
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D

import numpy as np # TODO: Remove numpy dependency.

class ModelNode(ABC):
    """
    A node in a scene.
    """
    @abstractmethod
    def __init__(self, fill_color: Color, fill_style: FillStyleType) -> None:
        """
        Initialize a node.
        
        :param Color fill_color: The color of the node.
        :param FillStyleType fill_style: The fill style of the node.
        """
        self.fill_color = fill_color
        self.fill_style = fill_style

    @abstractmethod
    def visible(self) -> float:
        """
        Calculate the visibility of the node.
        
        :return: The visibility of the node.
        :rtype: float
        """
        ...
    
    @abstractmethod
    def position_for_sorting(self) -> np.ndarray:
        """
        Calculate the position of the node for sorting.
        
        :return: The position of the node for sorting.
        :rtype: np.ndarray
        """
        ...

    @abstractmethod
    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the node.
        
        :param int res: The resolution of the node.
        :return: The points on the surface of the node, with shape (res, 3).
        :rtype: np.ndarray
        """
        ...

    @abstractmethod 
    def point_is_inside(self, point: np.ndarray) -> bool:
        """
        Check whether a point is inside the node.
        
        :param np.ndarray point: The point to check, with shape (3,).
        :return: Whether the point is inside the node.
        :rtype: bool
        """
        ...

    @abstractmethod
    def intersects_with(self, other: "ModelNode") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ModelNode other: The other node to check.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        ...

class ModelSphere(ModelNode):
    """
    A node in a scene.
    """
    def __init__(self, geometry: Sphere, fill_color: Color, fill_style: FillStyleType) -> None:
        """
        Initialize a node.
        
        :param Sphere geometry: The geometry of the node.
        :param Color fill_color: The color of the node.
        :param FillStyleType fill_style: The fill style of the node.
        """
        super().__init__(fill_color, fill_style)
        self.geometry = geometry

    def visible(self) -> float:
        """
        Calculate the visibility of the node.
        
        :return: The visibility of the node.
        :rtype: float
        """
        return self.geometry.radius > 0.0
    
    def position_for_sorting(self) -> np.ndarray:
        """
        Calculate the position of the node for sorting.
        
        :return: The position of the node for sorting.
        :rtype: np.ndarray
        """
        return self.geometry.center
    
    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the node.
        
        :param int res: The resolution of the node.
        :return: The points on the surface of the node, with shape (res, 3).
        :rtype: np.ndarray
        """
        return self.geometry.generate_points_on_surface(res)
    
    def point_is_inside(self, point: np.ndarray) -> bool:
        """
        Check whether a point is inside the node.
        
        :param np.ndarray point: The point to check, with shape (3,).
        :return: Whether the point is inside the node.
        :rtype: bool
        """
        return self.geometry.point_is_inside(point)
    
    def intersects_with(self, other: "ModelNode") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ModelNode other: The other node to check.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        if isinstance(other, ModelSphere):
            return np.linalg.norm(self.geometry.center - other.geometry.center) <= self.geometry.radius + other.geometry.radius
        elif isinstance(other, ModelCylinder):
            line = Line3D(other.geometry.start, other.geometry.end) 
            return line.distance_to_line(self.geometry.center) <= self.geometry.radius
        else:
            raise ValueError(f"Unknown node type '{type(other)}'")

class ModelCylinder(ModelNode):
    """
    An edge between two nodes.
    """
    def __init__(self, geometry: Cylinder, fill_color: Color, fill_style: FillStyleType) -> None:
        """
        Initialize an edge.
        
        :param Cylinder geometry: The geometry of the edge.
        :param Color fill_color: The color of the edge.
        :param FillStyleType fill_style: The fill style of the edge.
        """
        super().__init__(fill_color, fill_style)
        self.geometry = geometry

    def visible(self) -> float:
        """
        Calculate the visibility of the node.
        
        :return: The visibility of the node.
        :rtype: float
        """
        return self.geometry.radius > 0.0
    
    def position_for_sorting(self) -> np.ndarray:
        """
        Calculate the position of the node for sorting.
        
        :return: The position of the node for sorting.
        :rtype: np.ndarray
        """
        return (self.geometry.start + self.geometry.end) / 2
    
    def generate_points_on_surface(self, res: int) -> np.ndarray:
        """
        Generate points on the surface of the node.
        
        :param int res: The resolution of the node.
        :return: The points on the surface of the node, with shape (res, 3).
        :rtype: np.ndarray
        """
        return self.geometry.generate_points_on_surface(res)
    
    def point_is_inside(self, point: np.ndarray) -> bool:
        """
        Check whether a point is inside the node.
        
        :param np.ndarray point: The point to check, with shape (3,).
        :return: Whether the point is inside the node.
        :rtype: bool
        """
        return self.geometry.point_is_inside(point)
    
    def intersects_with(self, other: "ModelNode") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ModelNode other: The other node to check.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        line = Line3D(self.geometry.start, self.geometry.end) 
        if isinstance(other, ModelSphere):
            return line.distance_to_line(other.geometry.center) <= other.geometry.radius
        elif isinstance(other, ModelCylinder):
            # TODO: This is now a hack to sort of get bonds that are in the neighborhood.
            other_middle = (other.geometry.start + other.geometry.end) / 2
            return line.distance_to_line(other_middle) <= other.geometry.radius * 10
        else:
            raise ValueError(f"Unknown node type '{type(other)}'")

def get_node_polygon_vertices(this: ModelNode, others: ty.List[ModelNode], resolution: int) -> np.ndarray:
    """
    Get the vertices of the polygon that represents the visible part of the node.
    
    :param ModelNode this: The node to get the vertices of.
    :param ty.List[ModelNode] others: The nodes that the node intersects with.
    :param int resolution: The resolution of the node.
    :return: The vertices of the polygon that represents the visible part of the node.
    :rtype: np.ndarray
    """
    def is_visible(point: np.ndarray) -> bool:
        return all([not node.point_is_inside(point) for node in others])

    points = this.generate_points_on_surface(resolution) # [N, 3]
    visible = np.array([point[:2] for point in points if is_visible(point)])

    if len(visible) > 0:
        # inds = quick_hull_2d(visible)
        # verts = [visible[ind] for ind in inds]

        visible = [Point2D(point[0], point[1]) for point in visible]
        inds = calculate_convex_hull(visible)
        verts = [visible[ind] for ind in inds]
        verts = np.array([[vert.x, vert.y] for vert in verts])
        
        return verts
    else:
        return []

class Scene:
    def __init__(self, nodes: ty.List[ModelNode] = []) -> None:
        """
        Initialize a scene.

        :param ty.List[ModelNode] nodes: The nodes of the scene.
        """
        if not all(isinstance(node, ModelNode) for node in nodes):
            raise ValueError("All scene nodes must be of type ModelNode")

        self.nodes = nodes

    def __str__(self) -> str:
        """
        Return a string representation of the scene.
        
        :return: The string representation of the scene.
        :rtype: str
        """
        return f"Scene(nodes={len(self.nodes)})"

    def add_node(self, node: ModelNode) -> None:
        """
        Add a node to the scene.
        
        :param ModelNode node: The node to add.
        :return: None
        :rtype: None
        """
        self.nodes.append(node)

    def calculate_view_box(self, points: np.ndarray, margin: float = 5) -> ViewBox:
        """
        Calculate the view box of the scene.
        
        :param np.ndarray points: The points of the scene, with shape (num_points, 2).
        :param float margin: The margin around the scene.
        :return: The view box of the scene.
        :rtype: ViewBox
        """
        if len(points) == 0:
            ViewBox(0, 0, 0, 0)
        
        if points.shape != (len(points), 2):
            raise ValueError("Points must have shape (num_points, 2)")

        min_x = np.min(points[:, 0])
        min_y = np.min(points[:, 1])
        max_x = np.max(points[:, 0])
        max_y = np.max(points[:, 1])
    
        min_x -= margin
        min_y -= margin
        max_x += margin
        max_y += margin

        width = max_x - min_x
        height = max_y - min_y
        return ViewBox(min_x, min_y, width, height)

    def draw(self, res: int, verb: bool = False) -> str:
        """
        Draw the scene.

        :param int res: The resolution of the scene.
        :param bool verb: Whether to print progress.
        :return: The SVG string.
        :rtype: str
        """
        if not all(isinstance(node, ModelNode) for node in self.nodes):
            raise ValueError("All scene nodes must be of type ModelNode")
        
        nodes = [node for node in self.nodes if node.visible()]

        view_box = self.calculate_view_box(np.array([node.position_for_sorting()[:2] for node in nodes]))
        svg = Svg(view_box, version=1.0, encoding="UTF-8")

        nodes = sorted(nodes, key=lambda node: node.position_for_sorting()[2], reverse=False)

        objects, fills = [], []
        for i, node in enumerate(nodes):
            previous_nodes = [previous_node for previous_node in nodes[:i] if node.intersects_with(previous_node)]

            # Filter previous nodes to see which ones intersect with the current node.
            # TODO 

            # Create outline.
            if isinstance(node, ModelSphere) and len(previous_nodes) == 0:
                outline = Circle2D(f"node-{i}", node.geometry.center[:2], node.geometry.radius)
                objects.append(outline)

            else:
                points = get_node_polygon_vertices(node, previous_nodes, res)
                polygon = Polygon2D(f"node-{i}", points)
                objects.append(polygon)

            # Create outline fill.
            if node.fill_style == FillStyleType.Cartoon:
                fill = Fill(f"node-{i}", node.fill_color, Solid())
                fills.append(fill)

            elif node.fill_style == FillStyleType.Glossy:
                
                if isinstance(node, ModelSphere):
                    gradient = RadialGradient(node.geometry.center[:2], node.geometry.radius)
                    fill = Fill(f"node-{i}", node.fill_color, gradient)
                    fills.append(fill)

                elif isinstance(node, ModelCylinder):
                    gradient = LinearGradient(node.geometry.start[:2], node.geometry.end[:2], node.geometry.radius)
                    fill = Fill(f"node-{i}", node.fill_color, gradient)
                    fills.append(fill)

                else:
                    raise ValueError(f"Unknown node type '{type(node)}'")

            else:
                raise ValueError(f"Unknown fill style '{node.fill_style}'")

            # Only print if verbose is set to True.
            if verb:
                padding = len(str(len(nodes)))
                print(f"{i}".zfill(padding), end="\r")

        svg_str = svg.to_svg(fills, objects)

        return svg_str