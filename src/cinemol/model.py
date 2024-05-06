# -*- coding: utf-8 -*-

"""This model module contains classes for creating a scene with nodes and edges."""

import typing as ty
from logging import getLogger

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Cylinder,
    Line3D,
    Point2D,
    Point3D,
    Sphere,
    cylinder_intersects_with_cylinder,
    get_points_on_surface_cylinder,
    get_points_on_surface_sphere,
    point_is_inside_cylinder,
    point_is_inside_sphere,
    sphere_intersects_with_cylinder,
    sphere_intersects_with_sphere,
)
from cinemol.style import (
    Cartoon,
    Color,
    Depiction,
    Fill,
    FillStyle,
    Glossy,
    LinearGradient,
    RadialGradient,
    Solid,
    Wire,
)
from cinemol.svg import Line2D, Polygon2D, Shape2D, Svg, ViewBox

# ==============================================================================
# Model nodes
# ==============================================================================


class ModelSphere:
    """Sphere shaped node in a scene.

    Used as node in a scene.
    """

    def __init__(self, geometry: Sphere, depiction: Depiction) -> None:
        """Create a sphere shaped node in a scene.

        :param geometry: The geometry of the node.
        :type geometry: Sphere
        :param depiction: The depiction of the node.
        :type depiction: Depiction
        """
        self.geometry = geometry
        self.depiction = depiction


class ModelCylinder:
    """Cylinder shaped node in a scene.

    Used as edge between two nodes.
    """

    def __init__(self, geometry: Cylinder, depiction: Depiction) -> None:
        """Create a cylinder shaped node in a scene.

        :param geometry: The geometry of the edge.
        :type geometry: Cylinder
        :param depiction: The depiction of the edge.
        :type depiction: Depiction
        """
        self.geometry = geometry
        self.depiction = depiction


class ModelWire:
    """Wire shaped node in a scene.

    Used as edge between two nodes.
    """

    def __init__(self, geometry: Line3D, color: Color, width: float, opacity: float) -> None:
        """Create a wire shaped node in a scene.

        :param geometry: The geometry of the edge.
        :type geometry: Line3D
        :param color: The color of the edge.
        :type color: Color
        :param width: The width of the edge.
        :type width: float
        :param opacity: The opacity of the edge.
        :type opacity: float
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
    resolution: int,
    focal_length: ty.Optional[float] = None,
) -> ty.List[Point2D]:
    """Get the vertices of the polygon that represents the visible part of the node.

    :param this: The node to get the vertices of.
    :type this: ty.Union[ModelSphere, ModelCylinder, ModelWire]
    :param others: The other nodes in the scene.
    :type others: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
    :param resolution: The resolution of the polygon.
    :type resolution: int
    :param focal_length: The scaling factor of the polygon.
    :type focal_length: ty.Optional[float]
    :return: The vertices of the polygon.
    :rtype: ty.List[Point2D]
    """
    # Generate points on surface of node geometry.
    if isinstance(this, ModelSphere):
        points = get_points_on_surface_sphere(
            this.geometry, resolution, resolution, filter_for_pov=True
        )

    elif isinstance(this, ModelCylinder):
        points = get_points_on_surface_cylinder(this.geometry, int(resolution // 2.0))

    else:
        # If node is not a sphere or cylinder (i.e., unsupported geometries), return empty list.
        return []

    # Check if point is visible (i.e, not inside any other node geometry).
    visible_points = []
    for point in points:
        for node in others:
            if isinstance(node, ModelSphere) and point_is_inside_sphere(node.geometry, point):
                break

            elif isinstance(node, ModelCylinder) and point_is_inside_cylinder(node.geometry, point):
                break

        else:
            x, y, z = point.x, point.y, point.z

            if focal_length is not None:
                factor = focal_length / (z - focal_length)
                if factor < 0:  # Point is behind the point of view.
                    continue
            else:
                factor = 1.0

            visible_points.append(Point2D(x * factor, y * factor))

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


def prepare_nodes_for_intersecting(
    nodes_to_sort: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]],
    include_spheres: bool,
    include_cylinders: bool,
    include_wires: bool,
    rotation_over_x_axis: float = 0.0,
    rotation_over_y_axis: float = 0.0,
    rotation_over_z_axis: float = 0.0,
    scale: ty.Optional[float] = None,
) -> ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]:
    """Filter, rotate, and scale nodes based on what to include in the scene.

    :param nodes_to_sort: The nodes to filter.
    :type nodes_to_sort: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
    :param include_spheres: Whether to include spheres in the scene.
    :type include_spheres: bool
    :param include_cylinders: Whether to include cylinders in the scene.
    :type include_cylinders: bool
    :param include_wires: Whether to include wires in the scene.
    :type include_wires: bool
    :param rotation_over_x_axis: The rotation over the x-axis.
    :type rotation_over_x_axis: float
    :param rotation_over_y_axis: The rotation over the y-axis.
    :type rotation_over_y_axis: float
    :param rotation_over_z_axis: The rotation over the z-axis.
    :type rotation_over_z_axis: float
    :param scale: The scale of the scene.
    :type scale: ty.Optional[float]
    :return: The filtered nodes.
    :rtype: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
    """
    nodes: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]] = []

    for node in nodes_to_sort:
        if isinstance(node, ModelSphere) and include_spheres:
            node.geometry.center = node.geometry.center.rotate(
                rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis
            )

            if scale is not None:
                node.geometry.radius *= scale
                node.geometry.center = Point3D(
                    node.geometry.center.x * scale,
                    node.geometry.center.y * scale,
                    node.geometry.center.z * scale,
                )

                if isinstance(node.depiction, Cartoon):
                    node.depiction.outline_width *= scale

            nodes.append(node)

        elif isinstance(node, ModelCylinder) and include_cylinders:
            node.geometry.start = node.geometry.start.rotate(
                rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis
            )
            node.geometry.end = node.geometry.end.rotate(
                rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis
            )

            if scale is not None:
                node.geometry.radius *= scale
                node.geometry.start = Point3D(
                    node.geometry.start.x * scale,
                    node.geometry.start.y * scale,
                    node.geometry.start.z * scale,
                )
                node.geometry.end = Point3D(
                    node.geometry.end.x * scale,
                    node.geometry.end.y * scale,
                    node.geometry.end.z * scale,
                )

                if isinstance(node.depiction, Cartoon):
                    node.depiction.outline_width *= scale

            nodes.append(node)

        elif isinstance(node, ModelWire) and include_wires:
            node.geometry.start = node.geometry.start.rotate(
                rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis
            )
            node.geometry.end = node.geometry.end.rotate(
                rotation_over_x_axis, rotation_over_y_axis, rotation_over_z_axis
            )

            if scale is not None:
                node.geometry.start = Point3D(
                    node.geometry.start.x * scale,
                    node.geometry.start.y * scale,
                    node.geometry.start.z * scale,
                )
                node.geometry.end = Point3D(
                    node.geometry.end.x * scale,
                    node.geometry.end.y * scale,
                    node.geometry.end.z * scale,
                )

                node.width *= scale

            nodes.append(node)

    return nodes


def create_fill(
    node: ty.Union[ModelSphere, ModelCylinder, ModelWire], reference: str
) -> ty.Optional[Fill]:
    """Create a fill for a node.

    :param node: The node to create a fill for.
    :type node: ty.Union[ModelSphere, ModelCylinder, ModelWire]
    :param reference: The reference of the fill.
    :type reference: str
    :return: The fill.
    :rtype: ty.Optional[Fill]
    :raises AssertionError: If the geometry of the node is not appropriate for the
        node type.
    """
    fill = None

    if isinstance(node, ModelWire):
        stroke_color = node.color
        stroke_width = node.width
        opacity = node.opacity
        style: FillStyle = Wire(stroke_color, stroke_width, opacity)
        fill = Fill(reference, style)

    elif isinstance(node.depiction, Cartoon):
        fill_color = node.depiction.fill_color
        stroke_color = node.depiction.outline_color
        stroke_width = node.depiction.outline_width
        opacity = node.depiction.opacity
        style = Solid(fill_color, stroke_color, stroke_width, opacity)
        fill = Fill(reference, style)

    # Glossy style is different for spheres and cylinders.
    elif (
        isinstance(node.depiction, Glossy)
        and isinstance(node, ModelSphere)
        and isinstance(node.geometry, Sphere)
    ):
        if not isinstance(node.geometry, Sphere):
            raise AssertionError("Node geometry of ModelSphere must be a sphere.")

        x, y = node.geometry.center.x, node.geometry.center.y
        fill_color = node.depiction.fill_color
        center = Point2D(x, y)
        radius = node.geometry.radius
        opacity = node.depiction.opacity
        style = RadialGradient(fill_color, center, radius, opacity)
        fill = Fill(reference, style)

    elif isinstance(node.depiction, Glossy) and isinstance(node, ModelCylinder):
        if not isinstance(node.geometry, Cylinder):
            raise AssertionError("Node geometry of ModelCylinder must be a cylinder.")

        start_x, start_y = node.geometry.start.x, node.geometry.start.y
        end_x, end_y = node.geometry.end.x, node.geometry.end.y
        fill_color = node.depiction.fill_color
        start_center = Point2D(start_x, start_y)
        end_center = Point2D(end_x, end_y)
        radius = node.geometry.radius
        opacity = node.depiction.opacity
        style = LinearGradient(fill_color, start_center, end_center, opacity)
        fill = Fill(reference, style)

    return fill


def calculate_intersecting_nodes(
    node: ty.Union[ModelSphere, ModelCylinder, ModelWire],
    other_nodes: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]],
    calculate_sphere_sphere_intersections: bool,
    calculate_sphere_cylinder_intersections: bool,
    calculate_cylinder_sphere_intersections: bool,
    calculate_cylinder_cylinder_intersections: bool,
    filter_nodes_for_intersecting: bool,
) -> ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]:
    """Calculate which of the previous nodes intersect with the current node.

    :param node: The current node.
    :type node: ty.Union[ModelSphere, ModelCylinder, ModelWire]
    :param other_nodes: The previous nodes.
    :type other_nodes: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
    :param calculate_sphere_sphere_intersections: Whether to calculate intersections
        between spheres.
    :type calculate_sphere_sphere_intersections: bool
    :param calculate_sphere_cylinder_intersections: Whether to calculate intersections
        between spheres and cylinders.
    :type calculate_sphere_cylinder_intersections: bool
    :param calculate_cylinder_sphere_intersections: Whether to calculate intersections
        between cylinders and spheres.
    :type calculate_cylinder_sphere_intersections: bool
    :param calculate_cylinder_cylinder_intersections: Whether to calculate intersections
        between cylinders.
    :type calculate_cylinder_cylinder_intersections: bool
    :param filter_nodes_for_intersecting: Whether to filter nodes for intersecting nodes.
    :type filter_nodes_for_intersecting: bool
    :return: The previous nodes that intersect with the current node.
    :rtype: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
    """
    previous_nodes = []

    # Wireframe is drawn as a line and has no intersections with other nodes.
    if not isinstance(node, ModelWire):
        for prev_node in other_nodes:

            if not filter_nodes_for_intersecting:
                previous_nodes.append(prev_node)
                continue

            if (
                isinstance(node, ModelSphere)
                and isinstance(prev_node, ModelSphere)
                and calculate_sphere_sphere_intersections
            ):
                if sphere_intersects_with_sphere(node.geometry, prev_node.geometry):
                    previous_nodes.append(prev_node)

            elif (
                isinstance(node, ModelSphere) and isinstance(prev_node, ModelCylinder)
            ) and calculate_sphere_cylinder_intersections:
                if sphere_intersects_with_cylinder(node.geometry, prev_node.geometry):
                    previous_nodes.append(prev_node)

            elif (
                isinstance(node, ModelCylinder) and isinstance(prev_node, ModelSphere)
            ) and calculate_cylinder_sphere_intersections:
                if sphere_intersects_with_cylinder(prev_node.geometry, node.geometry):
                    previous_nodes.append(prev_node)

            elif (
                isinstance(node, ModelCylinder)
                and isinstance(prev_node, ModelCylinder)
                and calculate_cylinder_cylinder_intersections
            ):
                if cylinder_intersects_with_cylinder(node.geometry, prev_node.geometry):
                    previous_nodes.append(prev_node)

    return previous_nodes


class Scene:
    """A scene contains a list of nodes."""

    def __init__(
        self, nodes: ty.Optional[ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]] = None
    ) -> None:
        """Create a scene with nodes.

        :param nodes: The nodes in the scene.
        :type nodes: ty.List[ty.Union[ModelSphere, ModelCylinder, ModelWire]]
        """
        if nodes is None:
            nodes = []

        self.nodes = nodes

    def __str__(self) -> str:
        """Return a string representation of the scene.

        :return: The string representation of the scene.
        :rtype: str
        """
        return f"Scene(nodes={len(self.nodes)})"

    def add_node(self, node: ty.Union[ModelSphere, ModelCylinder, ModelWire]) -> None:
        """Add a node to the scene.

        :param node: The node to add.
        :type node: ty.Union[ModelSphere, ModelCylinder, ModelWire]
        """
        self.nodes.append(node)

    def calculate_view_box(self, points: ty.List[Point2D], margin: float) -> ViewBox:
        """Calculate the view box of the scene.

        :param points: The points to calculate the view box of.
        :type points: ty.List[Point2D]
        :param margin: The margin around the scene.
        :type margin: float
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
        focal_length: ty.Optional[float] = None,
    ) -> Svg:
        """
        Draw the scene.

        :param resolution: The resolution of the scene.
        :type resolution: int
        :param window: The window of the scene.
        :type window: ty.Optional[ty.Tuple[float, float]]
        :param view_box: The view box of the scene.
        :type view_box: ty.Optional[ViewBox]
        :param rotation_over_x_axis: The rotation over the x-axis.
        :type rotation_over_x_axis: float
        :param rotation_over_y_axis: The rotation over the y-axis.
        :type rotation_over_y_axis: float
        :param rotation_over_z_axis: The rotation over the z-axis.
        :type rotation_over_z_axis: float
        :param include_spheres: Whether to include spheres in the scene.
        :type include_spheres: bool
        :param include_cylinders: Whether to include cylinders in the scene.
        :type include_cylinders: bool
        :param include_wires: Whether to include wires in the scene.
        :type include_wires: bool
        :param calculate_sphere_sphere_intersections: Whether to calculate intersections
            between spheres.
        :type calculate_sphere_sphere_intersections: bool
        :param calculate_sphere_cylinder_intersections: Whether to calculate intersections
            between spheres and cylinders.
        :type calculate_sphere_cylinder_intersections: bool
        :param calculate_cylinder_sphere_intersections: Whether to calculate intersections
            between cylinders and spheres.
        :type calculate_cylinder_sphere_intersections: bool
        :param calculate_cylinder_cylinder_intersections: Whether to calculate intersections
            between cylinders.
        :type calculate_cylinder_cylinder_intersections: bool
        :param filter_nodes_for_intersecting: Whether to filter nodes for intersecting nodes.
        :type filter_nodes_for_intersecting: bool
        :param scale: The scale of the scene.
        :type scale: float
        :param focal_length: The focal length of the depiction. If None, the focal length
            is calculated based on the dimensions of the scene.
        :type focal_length: ty.Optional[float]
        :return: The scene as Svg.
        :rtype: Svg
        """
        logger = getLogger(__name__)

        # Filter geometries.
        nodes = prepare_nodes_for_intersecting(
            self.nodes,
            include_spheres,
            include_cylinders,
            include_wires,
            rotation_over_x_axis,
            rotation_over_y_axis,
            rotation_over_z_axis,
            scale,
        )

        # Get sorting values for nodes. We sort on z-coordinate as we always look at the
        # scene from the z-axis, towards the origin.
        sorting_values = []
        for node in nodes:
            if isinstance(node, ModelSphere):
                sorting_values.append(node.geometry.center.z)

            elif isinstance(node, ModelCylinder):
                cylinder_start = node.geometry.start
                cylinder_end = node.geometry.end
                midpoint_z = (cylinder_start.z + cylinder_end.z) / 2
                sorting_values.append(midpoint_z)

            elif isinstance(node, ModelWire):
                wire_start = node.geometry.start
                wire_end = node.geometry.end
                midpoint_z = (wire_start.z + wire_end.z) / 2
                sorting_values.append(midpoint_z)

        # Sort nodes by sorting values.
        nodes = [
            node
            for _, node in sorted(zip(sorting_values, nodes), key=lambda x: x[0], reverse=False)
        ]

        # Keep track of reference points for determining viewbox later on.
        if view_box is None:
            ref_points = []

        # Calculate 2D shape and fill for each node.
        objects: ty.List[Shape2D] = []
        fills = []

        # Loop over nodes and create shapes and fills to populate the SVG.
        for i, node in enumerate(nodes):

            # Create reference tag for node to connect shape to style.
            reference = f"node-{i}"

            # Calculate which of the previously drawn nodes intersect with the current node.
            previous_nodes = calculate_intersecting_nodes(
                node,
                nodes[:i],
                calculate_sphere_sphere_intersections,
                calculate_sphere_cylinder_intersections,
                calculate_cylinder_sphere_intersections,
                calculate_cylinder_cylinder_intersections,
                filter_nodes_for_intersecting,
            )

            # Calculate line for wire.
            if isinstance(node, ModelWire):
                start_x, start_y = node.geometry.start.x, node.geometry.start.y
                end_x, end_y = node.geometry.end.x, node.geometry.end.y

                if focal_length is not None:
                    start_z, end_z = node.geometry.start.z, node.geometry.end.z
                    start_factor = focal_length / (start_z - focal_length)
                    end_factor = focal_length / (end_z - focal_length)
                    start_x *= start_factor
                    start_y *= start_factor
                    end_x *= end_factor
                    end_y *= end_factor

                start: Point2D = Point2D(start_x, start_y)
                end: Point2D = Point2D(end_x, end_y)
                line = Line2D(reference, start, end)

                if view_box is None:
                    ref_points.extend([start, end])

                objects.append(line)

            # Otherwise, calculate polygon for visible part of node.
            else:
                points = get_node_polygon_vertices(node, previous_nodes, resolution, focal_length)
                polygon = Polygon2D(reference, points)

                if view_box is None:
                    ref_points.extend(points)

                objects.append(polygon)

            # Create style.
            fill = create_fill(node, reference)
            if fill is not None:
                fills.append(fill)

            padding = len(str(len(nodes)))
            logger.info(f" Drawing node {i + 1:>{padding}} of {len(nodes)}")

        # Calculate view box.
        if view_box is None:
            view_box = self.calculate_view_box(ref_points, 5.0)

        svg = Svg(view_box=view_box, window=window, fills=fills, objects=objects)

        return svg
