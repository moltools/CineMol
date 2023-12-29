"""
Contains definitions for a model scene.
"""
from abc import ABC, abstractmethod
import typing as ty

from cinemol.geometry import Point2D, Point3D
from cinemol.fitting import calculate_convex_hull
from cinemol.shapes import Line3D, Sphere, Cylinder 
from cinemol.style import Color, Fill, FillStyleType, Solid, RadialGradient, LinearGradient
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D

class ModelNode(ABC):
    """
    A node in a scene.
    """
    @abstractmethod
    def __init__(self, fill_color: Color, fill_style: FillStyleType) -> None:
        self.fill_color = fill_color
        self.fill_style = fill_style

    @abstractmethod
    def visible(self) -> float:
        ...
    
    @abstractmethod
    def position_for_sorting(self):
        ...

    @abstractmethod
    def generate_points_on_surface(self, res: int):
        ...

    @abstractmethod 
    def point_is_inside(self, point) -> bool:
        ...

    @abstractmethod
    def intersects_with(self, other: "ModelNode") -> bool:
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
    
    def position_for_sorting(self):
        """
        Calculate the position of the node for sorting.
        
        :return: The position of the node for sorting.
        :rtype: Point3D
        """
        return self.geometry.center
    
    def generate_points_on_surface(self, res: int):
        """
        Generate points on the surface of the node.
        
        :param int res: The resolution of the node.
        :return: The points on the surface of the node, with shape (res, 3).
        :rtype: ty.List[Point3D]
        """
        return self.geometry.generate_points_on_surface(res)
    
    def point_is_inside(self, point) -> bool:
        """
        Check whether a point is inside the node.
        
        :param Point3D point: The point to check, with shape (3,).
        :return: Whether the point is inside the node.
        :rtype: bool
        """
        if not isinstance(point, Point3D):
            point = Point3D(point[0], point[1], point[2])
        return self.geometry.point_is_inside(point)
    
    def intersects_with(self, other: "ModelNode") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ModelNode other: The other node to check.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        if isinstance(other, ModelSphere):
            dist = self.geometry.center.calculate_distance(other.geometry.center)
            return dist <= self.geometry.radius + other.geometry.radius
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
    
    def position_for_sorting(self):
        """
        Calculate the position of the node for sorting.
        
        :return: The position of the node for sorting.
        :rtype: Point3D
        """
        return Point3D(
            (self.geometry.start.x + self.geometry.end.x) / 2,
            (self.geometry.start.y + self.geometry.end.y) / 2,
            (self.geometry.start.z + self.geometry.end.z) / 2
        )
    
    def generate_points_on_surface(self, res: int):
        """
        Generate points on the surface of the node.
        
        :param int res: The resolution of the node.
        :return: The points on the surface of the node, with shape (res, 3).
        :rtype: ty.List[Point3D]
        """
        return self.geometry.generate_points_on_surface(res)
    
    def point_is_inside(self, point) -> bool:
        """
        Check whether a point is inside the node.
        
        :param Point3D point: The point to check, with shape (3,).
        :return: Whether the point is inside the node.
        :rtype: bool
        """
        if not isinstance(point, Point3D):
            point = Point3D(point[0], point[1], point[2])
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
            other_middle = other.geometry.start.midpoint(other.geometry.end)
            return line.distance_to_line(other_middle) <= other.geometry.radius * 10
        else:
            raise ValueError(f"Unknown node type '{type(other)}'")

def get_node_polygon_vertices(this: ModelNode, others: ty.List[ModelNode], resolution: int):
    """
    Get the vertices of the polygon that represents the visible part of the node.
    
    :param ModelNode this: The node to get the vertices of.
    :param ty.List[ModelNode] others: The nodes that the node intersects with.
    :param int resolution: The resolution of the node.
    :return: The vertices of the polygon that represents the visible part of the node.
    :rtype: ty.List[Point2D]
    """
    def is_visible(point) -> bool:
        return all([not node.point_is_inside(point) for node in others])

    points = this.generate_points_on_surface(resolution) # [N, 3]

    visible = [Point2D(point.x, point.y) for point in points if is_visible(point)]

    if len(visible) > 0:
        inds = calculate_convex_hull(visible)
        verts = [visible[ind] for ind in inds]
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

    def calculate_view_box(self, points: ty.List[Point2D], margin: float = 5) -> ViewBox:
        """
        Calculate the view box of the scene.
        
        :param ty.List[Point2D] points: The points to calculate the view box of.
        :param float margin: The margin around the scene.
        :return: The view box of the scene.
        :rtype: ViewBox
        """
        if len(points) == 0:
            ViewBox(0, 0, 0, 0)

        min_x, min_y, max_x, max_y = float("inf"), float("inf"), float("-inf"), float("-inf")
        for point in points:
            min_x = min(min_x, point.x)
            min_y = min(min_y, point.y)
            max_x = max(max_x, point.x)
            max_y = max(max_y, point.y)
    
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

        view_box = self.calculate_view_box([node.position_for_sorting() for node in nodes])
        svg = Svg(view_box, version=1.0, encoding="UTF-8")

        nodes = sorted(nodes, key=lambda node: node.position_for_sorting().z, reverse=False)

        has_spheres = any(isinstance(node, ModelSphere) for node in nodes)
        # if ther are spheres (even one), then don't calculate cylinder-cylinder intersections.

        objects, fills = [], []
        for i, node in enumerate(nodes):
            previous_nodes = [previous_node for previous_node in nodes[:i] if node.intersects_with(previous_node)]

            if isinstance(node, ModelSphere):
                previous_nodes = [node for node in previous_nodes if isinstance(node, ModelSphere)]
                # NOTE: otherwise some nodes disappear... spheres get cut off by cylinders that are behind it and intersect it.
            
            if has_spheres and isinstance(node, ModelCylinder):
                previous_nodes = [node for node in previous_nodes if isinstance(node, ModelSphere)]
                # NOTE: you can't see cylinders intersecting inside spheres anyway.

            # Filter previous nodes to see which ones intersect with the current node.
            # TODO 

            # Create outline.
            if isinstance(node, ModelSphere) and len(previous_nodes) == 0:
                outline = Circle2D(f"node-{i}", Point2D(node.geometry.center.x, node.geometry.center.y), node.geometry.radius)
                objects.append(outline)

            else:
                if isinstance(node, ModelCylinder):
                    temp_res = int(res / 4)
                else:
                    temp_res = res

                points = get_node_polygon_vertices(node, previous_nodes, temp_res)
                polygon = Polygon2D(f"node-{i}", points)
                objects.append(polygon)

            # Create outline fill.
            if node.fill_style == FillStyleType.Cartoon:
                fill = Fill(f"node-{i}", node.fill_color, Solid())
                fills.append(fill)

            elif node.fill_style == FillStyleType.Glossy:
                
                if isinstance(node, ModelSphere):
                    gradient = RadialGradient(
                        Point2D(node.geometry.center.x, node.geometry.center.y), 
                        node.geometry.radius
                    )
                    fill = Fill(f"node-{i}", node.fill_color, gradient)
                    fills.append(fill)

                elif isinstance(node, ModelCylinder):
                    gradient = LinearGradient(
                        Point2D(node.geometry.start.x, node.geometry.start.y),
                        Point2D(node.geometry.end.x, node.geometry.end.y), 
                        node.geometry.radius
                    )
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