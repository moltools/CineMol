"""
Contains definitions for a model scene.
"""
import typing as ty

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Point2D, 
    Point3D,
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
from cinemol.svg import ViewBox, Svg, Polygon2D, Line2D

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
    pov_z: float,
    others: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]], 
    resolution: int
) -> ty.List[Point2D]:
    """
    Get the vertices of the polygon that represents the visible part of the node.
    
    :param ty.Union[ModelSphere, ModelCylinder, ModelWire] this: The node to get the vertices of.
    :param float pov_z: The z-coordinate of the point of view.
    :param ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]] others: The other nodes in the scene.
    :param int resolution: The resolution of the polygon.
    :return: The vertices of the polygon.
    :rtype: ty.List[Point2D]
    """
    # Generate points on surface of node geometry.
    if isinstance(this, ModelSphere):
        points = get_points_on_surface_sphere(this.geometry, resolution, resolution, filter_for_pov=True)
    
    elif isinstance(this, ModelCylinder):
        points = get_points_on_surface_cylinder(this.geometry, int(resolution // 2.0))

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
            x, y = point.x, point.y
            visible_points.append(Point2D(x, y))

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
        window: ty.Optional[ty.Tuple[float, float]] = None,
        view_box: ty.Optional[ViewBox] = None,
        rotation_over_x_axis: float = 0.0,
        rotation_over_y_axis: float = 0.0,
        rotation_over_z_axis: float = 0.0,
        include_spheres: bool = True,
        include_cylinders: bool = True,
        include_wires: bool = True,
        calculate_sphere_sphere_intersections: bool = True,
        calculate_sphere_cylinder_intersections: bool = True,
        calculate_cylinder_sphere_intersections: bool = True,
        calculate_cylinder_cylinder_intersections: bool = True,
        filter_nodes_for_intersecting: bool = True,
        scale: float = 1.0,
        focal_length: float = 10.0,
        verbose: bool = False,
    ) -> Svg:
        """
        Draw the scene.

        :param int resolution: The resolution of the scene.
        :param ty.Optional[ty.Tuple[float, float]] window: The window of the scene.
        :param ty.Optional[ViewBox] view_box: The view box of the scene.
        :param float rotation_over_x_axis: The rotation over the x-axis.
        :param float rotation_over_y_axis: The rotation over the y-axis.
        :param float rotation_over_z_axis: The rotation over the z-axis.
        :param bool include_spheres: Whether to include spheres in the scene.
        :param bool include_cylinders: Whether to include cylinders in the scene.
        :param bool include_wires: Whether to include wires in the scene.
        :param bool calculate_sphere_sphere_intersections: Whether to calculate intersections between spheres.
        :param bool calculate_sphere_cylinder_intersections: Whether to calculate intersections between spheres and cylinders.
        :param bool calculate_cylinder_sphere_intersections: Whether to calculate intersections between cylinders and spheres.
        :param bool calculate_cylinder_cylinder_intersections: Whether to calculate intersections between cylinders.
        :param bool filter_nodes_for_intersecting: Whether to filter nodes for intersecting nodes.
        :param float scale: The scale of the scene.
        :param float focal_length: The focal length of the scene.
        :param bool verbose: Whether to print progress.
        :return: The scene as Svg.
        :rtype: Svg
        """
        # Filter geometries.
        nodes = []
        for node in self.nodes: 
            if isinstance(node, ModelSphere) and include_spheres:
                node.geometry.center = node.geometry.center.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                
                if scale is not None:
                    node.geometry.radius *= scale
                    node.geometry.center = Point3D(
                        node.geometry.center.x * scale,
                        node.geometry.center.y * scale,
                        node.geometry.center.z * scale
                    )

                    if isinstance(node.depiction, Cartoon):
                        node.depiction.outline_width *= scale

                nodes.append(node)
            
            elif isinstance(node, ModelCylinder) and include_cylinders:
                node.geometry.start = node.geometry.start.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                node.geometry.end = node.geometry.end.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)

                if scale is not None:
                    node.geometry.radius *= scale
                    node.geometry.start = Point3D(
                        node.geometry.start.x * scale,
                        node.geometry.start.y * scale,
                        node.geometry.start.z * scale
                    )
                    node.geometry.end = Point3D(
                        node.geometry.end.x * scale,
                        node.geometry.end.y * scale,
                        node.geometry.end.z * scale
                    )

                    if isinstance(node.depiction, Cartoon):
                        node.depiction.outline_width *= scale

                nodes.append(node)

            elif isinstance(node, ModelWire) and include_wires:
                node.geometry.start = node.geometry.start.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)
                node.geometry.end = node.geometry.end.rotate(rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis)

                if scale is not None:
                    node.geometry.start = Point3D(
                        node.geometry.start.x * scale,
                        node.geometry.start.y * scale,
                        node.geometry.start.z * scale
                    )
                    node.geometry.end = Point3D(
                        node.geometry.end.x * scale,
                        node.geometry.end.y * scale,
                        node.geometry.end.z * scale
                    )

                    node.width *= scale

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

        # Keep track of reference points for determining viewbox later on.
        if view_box is None:
            ref_points = []

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

                    if not filter_nodes_for_intersecting:
                        previous_nodes.append(prev_node)
                        continue

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

            # Calculate line for wire.
            if isinstance(node, ModelWire):
                start_x, start_y = node.geometry.start.x, node.geometry.start.y 
                end_x, end_y = node.geometry.end.x, node.geometry.end.y
                start = Point2D(start_x, start_y)
                end = Point2D(end_x, end_y)
                line = Line2D(reference, start, end)
                
                if view_box is None:
                    ref_points.extend([start, end])

                objects.append(line)
            
            # Otherwise, calculate polygon for visible part of node.
            else:
                pov_z = focal_length
                points = get_node_polygon_vertices(node, pov_z, previous_nodes, resolution)
                polygon = Polygon2D(reference, points)

                if view_box is None:
                    ref_points.extend(points)

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
                x, y = node.geometry.center.x, node.geometry.center.y
                fill_color = node.depiction.fill_color
                center = Point2D(x, y)
                radius = node.geometry.radius
                opacity = node.depiction.opacity
                style = RadialGradient(fill_color, center, radius, opacity)
                fills.append(Fill(reference, style))
            
            elif isinstance(node.depiction, Glossy) and isinstance(node, ModelCylinder):
                start_x, start_y = node.geometry.start.x, node.geometry.start.y
                end_x, end_y = node.geometry.end.x, node.geometry.end.y
                fill_color = node.depiction.fill_color
                start_center = Point2D(start_x, start_y)
                end_center = Point2D(end_x, end_y)
                radius = node.geometry.radius
                opacity = node.depiction.opacity
                style = LinearGradient(fill_color, start_center, end_center, opacity)
                fills.append(Fill(reference, style))
            
            # Only print if verbose is set to True.
            if verbose:
                padding = len(str(len(nodes)))
                print(f"{i}".zfill(padding), end="\r")

        # Calculate view box.
        if view_box is None:
            view_box = self.calculate_view_box(ref_points, 5.0)
            
        svg = Svg(view_box=view_box, window=window, fills=fills, objects=objects)

        return svg