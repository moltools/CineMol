"""
Contains definitions for a model scene.
"""
import typing as ty

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Point2D, 
    Sphere, 
    Cylinder,
    Line3D,
    sphere_intersects_with_sphere,
    sphere_intersects_with_cylinder,
    cylinder_intersects_with_cylinder,
    get_points_on_surface_sphere,
    get_points_on_surface_cylinder,
    point_is_inside_sphere,
    point_is_inside_cylinder
)
from cinemol.style import Color, Depiction, Cartoon, Glossy, Fill, Wire, Solid, RadialGradient, LinearGradient
from cinemol.svg import ViewBox, Svg, Circle2D, Polygon2D, Line2D

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

class ModelWire:
    """
    A wire between two nodes.
    """
    def __init__(self, geometry: Line3D, color: Color, width: float, opacity: float) -> None:
        """
        Initialize a wire.
        
        :param Line3D geometry: The geometry of the wire.
        :param Color color: The color of the wire.
        :param float width: The width of the wire.
        :param float opacity: The opacity of the wire.
        """
        self.geometry = geometry
        self.color = color
        self.width = width
        self.opacity = opacity

# ==============================================================================
# Create visible 2D polygon from node geometry
# ==============================================================================

def get_node_polygon_vertices(
    this: ty.Union[ModelSphere, ModelCylinder, ModelWire], 
    others: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]], 
    resolution: int
) -> ty.List[Point2D]:
    """
    Get the vertices of the polygon that represents the visible part of the node.
    
    :param ty.Union[ModelSphere, ModelCylinder, ModelWire] this: The node to get the vertices of.
    :param ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]] others: The other nodes in the scene.
    :param int resolution: The resolution of the polygon.
    :return: The vertices of the polygon.
    :rtype: ty.List[Point2D]
    """
    # Generate points on surface of node geometry.
    if isinstance(this, ModelSphere):
        points = get_points_on_surface_sphere(this.geometry, resolution, resolution, filter_for_pov=True)
    
    elif isinstance(this, ModelCylinder):
        points = get_points_on_surface_cylinder(this.geometry, resolution)

    else:
        # If node is not a sphere or cylinder (i.e., unsupported geometries), return empty list. 
        return []

    # Check if point is visible (i.e, not inside any other node geometry).
    visible_points = []
    for point in points:
        for node in others:
            if isinstance(node, ModelSphere) and point_is_inside_sphere(node.geometry, point): break 
            elif isinstance(node, ModelCylinder) and point_is_inside_cylinder(node.geometry, point): break
        else:
            # Convert point to 2D point by ignoring z-coordinate.
            visible_points.append(Point2D(point.x, point.y))    

    # If no visible points, return empty list.
    if len(visible_points) == 0: 
        return []

    # Calculate convex hull of visible points.
    inds = calculate_convex_hull(visible_points)
    verts = [visible_points[ind] for ind in inds]
    
    return verts

# ==============================================================================
# Draw scene
# ==============================================================================

class Scene:
    def __init__(self, nodes: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]] = []) -> None:
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

    def add_node(self, node: ty.Union[ModelSphere, ModelCylinder, ModelWire]) -> None:
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
        rotation_over_x_axis: bool = 0.0,
        rotation_over_y_axis: bool = 0.0,
        rotation_over_z_axis: bool = 0.0,
        include_spheres: bool = True,
        include_cylinders: bool = True,
        include_wires: bool = True,
        calculate_sphere_sphere_intersections: bool = True,
        calculate_sphere_cylinder_intersections: bool = True,
        calculate_cylinder_sphere_intersections: bool = True,
        calculate_cylinder_cylinder_intersections: bool = True,
        svg_version: float = 1.0,
        svg_encoding: str = "UTF-8"
    ) -> str:
        """
        Draw the scene.

        :param int resolution: The resolution of the scene.
        :param bool verbose: Whether to print progress.
        :param float rotation_over_x_axis: The rotation over the x-axis.
        :param float rotation_over_y_axis: The rotation over the y-axis.
        :param float rotation_over_z_axis: The rotation over the z-axis.
        :param bool include_spheres: Whether to include spheres in the scene.
        :param bool include_cylinders: Whether to include cylinders in the scene.
        :param bool include_wires: Whether to include wires in the scene.
        :param bool calculate_sphere_spere_intersections: Whether to calculate sphere-sphere intersections.
        :param bool calculate_sphere_cylinder_intersections: Whether to calculate sphere-cylinder intersections.
        :param bool calculate_cylinder_cylinder_intersections: Whether to calculate cylinder-cylinder intersections.
        :return: The SVG string.
        :rtype: str
        """
        # Filter geometries.
        nodes = []
        for node in self.nodes: 
            if isinstance(node, ModelSphere) and include_spheres:
                node.geometry.center = node.geometry.center.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                nodes.append(node)
            
            elif isinstance(node, ModelCylinder) and include_cylinders:
                node.geometry.start = node.geometry.start.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                node.geometry.end = node.geometry.end.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                nodes.append(node)

            elif isinstance(node, ModelWire) and include_wires:
                node.geometry.start = node.geometry.start.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                node.geometry.end = node.geometry.end.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
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

            elif isinstance(node, ModelWire):
                start, end = node.geometry.start, node.geometry.end
                midpoint_z = (start.z + end.z) / 2
                sorting_values.append(midpoint_z)

        # Sort nodes by sorting values.
        nodes = [
            node for _, node 
            in sorted(
                zip(sorting_values, nodes), 
                key=lambda x: x[0], reverse=False
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

            elif isinstance(node, ModelWire):
                start, end = node.geometry.start, node.geometry.end
                midpoint = Point2D((start.x + end.x) / 2, (start.y + end.y) / 2)
                points.append(midpoint)

        view_box = self.calculate_view_box(points, margin=5)
        svg = Svg(view_box, svg_version, svg_encoding)

        # Calculate 2D shape and fill for each node.
        objects, fills = [], []

        for i, node in enumerate(nodes):
            
            # Create reference tag for node to connect shape to style.
            reference = f"node-{i}"

            # Calculate which of the previously drawn nodes intersect with the current node.
            previous_nodes = []

            # Wireframe is drawn as a line and has no intersections with previous nodes.
            if not isinstance(node, ModelWire):
                for prev_node in nodes[:i]:

                    if (
                        isinstance(node, ModelSphere) and 
                        isinstance(prev_node, ModelSphere) and 
                        calculate_sphere_sphere_intersections
                    ):
                        if sphere_intersects_with_sphere(node.geometry, prev_node.geometry):
                            previous_nodes.append(prev_node)

                    elif (
                        isinstance(node, ModelSphere) and 
                        isinstance(prev_node, ModelCylinder)
                    ) and calculate_sphere_cylinder_intersections:
                        if sphere_intersects_with_cylinder(node.geometry, prev_node.geometry):
                            previous_nodes.append(prev_node)
                
                    elif (
                        isinstance(node, ModelCylinder) and 
                        isinstance(prev_node, ModelSphere)
                    ) and calculate_cylinder_sphere_intersections:
                        if sphere_intersects_with_cylinder(prev_node.geometry, node.geometry):
                            previous_nodes.append(prev_node)

                    elif (
                        isinstance(node, ModelCylinder) and 
                        isinstance(prev_node, ModelCylinder) and 
                        calculate_cylinder_cylinder_intersections
                    ):
                        if cylinder_intersects_with_cylinder(node.geometry, prev_node.geometry):
                            previous_nodes.append(prev_node)     

            if isinstance(node, ModelWire):
                start = Point2D(node.geometry.start.x, node.geometry.start.y)
                end = Point2D(node.geometry.end.x, node.geometry.end.y)
                line = Line2D(reference, start, end)
                objects.append(line)

            # Create outline for model spheres with no intersections with previous nodes.
            elif isinstance(node, ModelSphere) and len(previous_nodes) == 0:
                center = Point2D(node.geometry.center.x, node.geometry.center.y)
                outline = Circle2D(reference, center, node.geometry.radius)
                objects.append(outline)
            
            # Otherwise, calculate polygon for visible part of node.
            else:
                points = get_node_polygon_vertices(node, previous_nodes, resolution)
                polygon = Polygon2D(reference, points)
                objects.append(polygon)

            # Create style.
            if isinstance(node, ModelWire):
                stroke_color = node.color
                stroke_width = node.width
                opacity = node.opacity
                style = Wire(stroke_color, stroke_width, opacity)
                fills.append(Fill(reference, style))

            elif isinstance(node.depiction, Cartoon):
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
                style = LinearGradient(fill_color, start_center, end_center, opacity)
                fills.append(Fill(reference, style))
            
            # Only print if verbose is set to True.
            if verbose:
                padding = len(str(len(nodes)))
                print(f"{i}".zfill(padding), end="\r")

        # Draw the scene.
        svg_str = svg.to_svg(fills, objects)

        return svg_str