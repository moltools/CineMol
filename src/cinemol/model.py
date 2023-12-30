"""
Contains definitions for a model scene.
"""
import typing as ty

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Point2D, 
    Line3D, 
    Sphere, 
    CylinderCapType,
    Cylinder,
    distance_to_line
)
from cinemol.style import Color, Fill, FillStyleType, Solid, RadialGradient, LinearGradient
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D

# ==============================================================================
# Model nodes
# ==============================================================================

class ModelSphere:
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
        self.geometry = geometry
        self.fill_color = fill_color
        self.fill_style = fill_style
    
    def intersects_with(self, other: "ty.Union[ModelSphere, ModelCylinder]") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ty.Union[ModelSphere, ModelCylinder] other: The other node.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        raise NotImplementedError("TODO")

class ModelCylinder:
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
        self.geometry = geometry
        self.fill_color = fill_color
        self.fill_style = fill_style
    
    def intersects_with(self, other: "ty.Union[ModelSphere, ModelCylinder]") -> bool:
        """
        Check whether the node intersects with another node.
        
        :param ty.Union[ModelSphere, ModelCylinder] other: The other node.
        :return: Whether the node intersects with the other node.
        :rtype: bool
        """
        raise NotImplementedError("TODO")

# ==============================================================================
# Create visible polygon
# ==============================================================================

def get_node_polygon_vertices(
    this: ty.Union[ModelSphere, ModelCylinder], 
    others: ty.List[ty.Union[ModelSphere, ModelCylinder]], 
    resolution: int
) -> ty.List[Point2D]:
    """
    Get the vertices of the polygon that represents the visible part of the node.
    
    :param ty.Union[ModelSphere, ModelCylinder] this: The node to get the vertices of.
    :param ty.List[ty.Union[ModelSphere, ModelCylinder]] others: The other nodes in the scene.
    :param int resolution: The resolution of the polygon.
    :return: The vertices of the polygon.
    :rtype: ty.List[Point2D]
    """
    raise NotImplementedError("TODO")

    # def is_visible(point) -> bool:
    #     return all([not node.point_is_inside(point) for node in others])

    # points = this.generate_points_on_surface(resolution) # [N, 3]

    # visible = [Point2D(point.x, point.y) for point in points if is_visible(point)]

    # if len(visible) > 0:
    #     inds = calculate_convex_hull(visible)
    #     verts = [visible[ind] for ind in inds]
    #     return verts
    # else:
    #     return []

# ==============================================================================
# Draw scene
# ==============================================================================

class Scene:
    def __init__(self, nodes: ty.List[ty.Union[ModelSphere, ModelCylinder]] = []) -> None:
        """
        Initialize a scene.

        :param ty.List[ty.Union[ModelSphere, ModelCylinder]] nodes: The nodes in the scene.
        """
        self.nodes = nodes

    def __str__(self) -> str:
        """
        Return a string representation of the scene.
        
        :return: The string representation of the scene.
        :rtype: str
        """
        return f"Scene(nodes={len(self.nodes)})"

    def add_node(self, node: ty.Union[ModelSphere, ModelCylinder]) -> None:
        """
        Add a node to the scene.
        
        :param ty.Union[ModelSphere, ModelCylinder] node: The node to add.
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

    def draw(self, resolution: int, verbose: bool = False) -> str:
        """
        Draw the scene.

        :param int resolution: The resolution of the scene.
        :param bool verbose: Whether to print progress.
        :return: The SVG string.
        :rtype: str
        """
        raise NotImplementedError("TODO")

        # if not all(isinstance(node, ModelNode) for node in self.nodes):
        #     raise ValueError("All scene nodes must be of type ModelNode")
        
        # nodes = [node for node in self.nodes if node.visible()]

        # view_box = self.calculate_view_box([node.position_for_sorting() for node in nodes])
        # svg = Svg(view_box, version=1.0, encoding="UTF-8")

        # nodes = sorted(nodes, key=lambda node: node.position_for_sorting().z, reverse=False)

        # has_spheres = any(isinstance(node, ModelSphere) for node in nodes)
        # # if ther are spheres (even one), then don't calculate cylinder-cylinder intersections.

        # objects, fills = [], []
        # for i, node in enumerate(nodes):
        #     previous_nodes = [previous_node for previous_node in nodes[:i] if node.intersects_with(previous_node)]

        #     if isinstance(node, ModelSphere):
        #         previous_nodes = [node for node in previous_nodes if isinstance(node, ModelSphere)]
        #         # NOTE: otherwise some nodes disappear... spheres get cut off by cylinders that are behind it and intersect it.
            
        #     if has_spheres and isinstance(node, ModelCylinder):
        #         previous_nodes = [node for node in previous_nodes if isinstance(node, ModelSphere)]
        #         # NOTE: you can't see cylinders intersecting inside spheres anyway.

        #     # Filter previous nodes to see which ones intersect with the current node.
        #     # TODO 

        #     # Create outline.
        #     if isinstance(node, ModelSphere) and len(previous_nodes) == 0:
        #         outline = Circle2D(f"node-{i}", Point2D(node.geometry.center.x, node.geometry.center.y), node.geometry.radius)
        #         objects.append(outline)

        #     else:
        #         if isinstance(node, ModelCylinder):
        #             temp_res = int(res / 4)
        #         else:
        #             temp_res = res

        #         points = get_node_polygon_vertices(node, previous_nodes, temp_res)
        #         polygon = Polygon2D(f"node-{i}", points)
        #         objects.append(polygon)

        #     # Create outline fill.
        #     if node.fill_style == FillStyleType.Cartoon:
        #         fill = Fill(f"node-{i}", node.fill_color, Solid())
        #         fills.append(fill)

        #     elif node.fill_style == FillStyleType.Glossy:
                
        #         if isinstance(node, ModelSphere):
        #             gradient = RadialGradient(
        #                 Point2D(node.geometry.center.x, node.geometry.center.y), 
        #                 node.geometry.radius
        #             )
        #             fill = Fill(f"node-{i}", node.fill_color, gradient)
        #             fills.append(fill)

        #         elif isinstance(node, ModelCylinder):
        #             gradient = LinearGradient(
        #                 Point2D(node.geometry.start.x, node.geometry.start.y),
        #                 Point2D(node.geometry.end.x, node.geometry.end.y), 
        #                 node.geometry.radius
        #             )
        #             fill = Fill(f"node-{i}", node.fill_color, gradient)
        #             fills.append(fill)

        #         else:
        #             raise ValueError(f"Unknown node type '{type(node)}'")

        #     else:
        #         raise ValueError(f"Unknown fill style '{node.fill_style}'")

        #     # Only print if verbose is set to True.
        #     if verb:
        #         padding = len(str(len(nodes)))
        #         print(f"{i}".zfill(padding), end="\r")

        # svg_str = svg.to_svg(fills, objects)

        # return svg_str