"""
Contains definitions for a model scene.
"""
import typing as ty

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Point2D, 
    Point3D,
    Line3D, 
    Sphere, 
    CylinderCapType,
    Cylinder,
    distance_to_line
)
from cinemol.style import Depiction, Cartoon, Glossy, Fill, Solid, RadialGradient, LinearGradient
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D

# ==============================================================================
# Model nodes
# ==============================================================================

class ModelSphere:
    """
    A node in a scene.
    """
    def __init__(self, geometry: Sphere, depiction: Depiction) -> None:
        """
        Initialize a node.
        
        :param Sphere geometry: The geometry of the node.
        :param Depiction depiction: The depiction of the node.
        """
        self.geometry = geometry
        self.depiction = depiction 
    
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
    def __init__(self, geometry: Cylinder, depiction: Depiction) -> None:
        """
        Initialize an edge.
        
        :param Cylinder geometry: The geometry of the edge.
        :param Depiction depiction: The depiction of the edge.
        """
        self.geometry = geometry
        self.depiction = depiction
    
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

    def calculate_view_box(self, points: ty.List[Point2D], margin: float) -> ViewBox:
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

    def draw(
        self, 
        resolution: int, 
        verbose: bool = False,
        include_spheres: bool = True,
        include_cylinders: bool = True,
        calculate_sphere_spere_intersections: bool = True,
        calculate_sphere_cylinder_intersections: bool = True,
        calculate_cylinder_cylinder_intersections: bool = True
    ) -> str:
        """
        Draw the scene.

        :param int resolution: The resolution of the scene.
        :param bool verbose: Whether to print progress.
        :param bool include_spheres: Whether to include spheres in the scene.
        :param bool include_cylinders: Whether to include cylinders in the scene.
        :param bool calculate_sphere_spere_intersections: Whether to calculate sphere-sphere intersections.
        :param bool calculate_sphere_cylinder_intersections: Whether to calculate sphere-cylinder intersections.
        :param bool calculate_cylinder_cylinder_intersections: Whether to calculate cylinder-cylinder intersections.
        :return: The SVG string.
        :rtype: str
        """
        # Make sure only sphere and cylinder geometries are in the scene.
        nodes = []
        for node in self.nodes: 
            if isinstance(node, ModelSphere) and include_spheres:
                nodes.append(node)
            
            elif isinstance(node, ModelCylinder) and include_cylinders:
                nodes.append(node)

        # Get sorting values for nodes. We sort on z-coordinate as we always look at the
        # scene from the z-axis, towards the origin.
        sorting_values = []
        for node in nodes:
            if isinstance(node, ModelSphere):
                sorting_values.append(node.geometry.center.z)

            elif isinstance(node, ModelCylinder):
                start, end = node.geometry.start, node.geometry.end
                midpoint_z = (start.z + end.z) / 2
                sorting_values.append(midpoint_z)

        # Sort nodes by sorting values.
        nodes = [
            node for _, node 
            in sorted(
                zip(sorting_values, nodes), 
                key=lambda x: x[0], reverse=True
            )
        ]

        # Calculate size of view box. Only use x and y coordinates of node geometries.
        points = []
        for node in nodes:
            if isinstance(node, ModelSphere):
                point = Point2D(node.geometry.center.x, node.geometry.center.y)
                points.append(point)
            
            elif isinstance(node, ModelCylinder):
                start, end = node.geometry.start, node.geometry.end
                midpoint = Point2D((start.x + end.x) / 2, (start.y + end.y) / 2)
                points.append(midpoint)

        view_box = self.calculate_view_box(points, margin=5)
        svg = Svg(view_box, version=1.0, encoding="UTF-8")

        # Calculate 2D shape and fill for each node.
        objects, fills = [], []

        for i, node in enumerate(nodes):
            
            # Create reference tag for node to connect shape to style.
            reference = f"node-{i}"

            # Calculate which of the previously drawn nodes intersect with the current node.
            previous_nodes = []
            for prev_node in nodes[:i]:

                if (
                    isinstance(node, ModelSphere) and 
                    isinstance(prev_node, ModelSphere) and 
                    calculate_sphere_spere_intersections
                ):
                    if node.intersects_with(prev_node):
                        previous_nodes.append(prev_node)

                elif (
                    isinstance(node, ModelSphere) and isinstance(prev_node, ModelCylinder) or
                    isinstance(node, ModelCylinder) and isinstance(prev_node, ModelSphere)
                ) and calculate_sphere_cylinder_intersections:
                    if node.intersects_with(prev_node):
                        previous_nodes.append(prev_node)

                elif (
                    isinstance(node, ModelCylinder) and 
                    isinstance(prev_node, ModelCylinder) and 
                    calculate_cylinder_cylinder_intersections
                ):
                    if node.intersects_with(prev_node):
                        previous_nodes.append(prev_node)

            # Create outline for model spheres with no intersections with previous nodes.
            if isinstance(node, ModelSphere) and len(previous_nodes) == 0:
                center = Point2D(node.geometry.center.x, node.geometry.center.y)
                outline = Circle2D(reference, center, node.geometry.radius)
                objects.append(outline)
            
            # Otherwise, calculate polygon for visible part of node.
            else:
                points = get_node_polygon_vertices(node, previous_nodes, resolution)
                polygon = Polygon2D(reference, points)
                objects.append(polygon)

            # Create fill for node.
            if isinstance(node.depiction, Cartoon):
                fill_color = node.depiction.fill_color
                stroke_color = node.depiction.outline_color
                stroke_width = node.depiction.outline_width
                opacity = node.depiction.opacity
                style = Solid(fill_color, stroke_color, stroke_width, opacity)
                fills.append(Fill(reference, style))

            # Glossy style is different for spheres and cylinders.
            elif isinstance(node.depiction, Glossy) and isinstance(node, ModelSphere):
                fill_color = node.depiction.fill_color
                center = Point2D(node.geometry.center.x, node.geometry.center.y)
                radius = node.geometry.radius
                opacity = node.depiction.opacity
                style = RadialGradient(fill_color, center, radius, opacity)
                fills.append(Fill(reference, style))
            
            elif isinstance(node.depiction, Glossy) and isinstance(node, ModelCylinder):
                fill_color = node.depiction.fill_color
                start_center = Point2D(node.geometry.start.x, node.geometry.start.y)
                end_center = Point2D(node.geometry.end.x, node.geometry.end.y)
                radius = node.geometry.radius
                opacity = node.depiction.opacity
                style = LinearGradient(fill_color, start_center, end_center, radius, opacity)
                fills.append(Fill(reference, style))
            
            # Only print if verbose is set to True.
            if verbose:
                padding = len(str(len(nodes)))
                print(f"{i}".zfill(padding), end="\r")

        # Draw the scene.
        svg_str = svg.to_svg(fills, objects)

        return svg_str