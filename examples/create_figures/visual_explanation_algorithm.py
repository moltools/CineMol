#!/usr/bin/env python3
"""
Description:    Draw step-wise explanation of the algorithm.
Dependencies:   matplotlib==3.8.4; numpy==1.26.4
Usage:          python3 visual_explanation_algorithm.py -o /path/to/output/dir [-d resolution] [--t]
"""
import argparse
import math
import os
import typing as ty

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from cinemol.fitting import calculate_convex_hull
from cinemol.geometry import (
    Circle3D,
    Cylinder,
    CylinderCapType,
    Line3D,
    Plane3D,
    Point2D,
    Point3D,
    Sphere,
    get_points_on_circumference_circle_3d,
    get_points_on_line_3d,
    point_is_inside_cylinder,
    point_is_inside_sphere,
    same_side_of_plane,
)

# ==============================================================================
# Command line interface
# ==============================================================================


def cli() -> argparse.Namespace:
    """
    Command line interface for this script.

    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=str, required=True, help="Path to output dir.")
    parser.add_argument("-d", type=int, required=False, default=300, help="DPI of output images.")
    parser.add_argument("--t", action="store_true", required=False, help="Transparent background.")
    return parser.parse_args()


# ==============================================================================
# Helper functions
# ==============================================================================


def get_axes_3d() -> Axes3D:
    """
    Get 3D axes.

    :return: 3D axes.
    :rtype: Axes3D
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Remove axes.
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # Plot 3D plot with three axes x,y,z as black arrows.
    line_length = 2.0
    line_color = "black"

    ax.plot([0, line_length], [0, 0], [0, 0], color=line_color, linewidth=2, zorder=1)
    ax.text(line_length + 0.1, -0.1, 0, "X", fontsize=12, color=line_color, zorder=1)

    ax.plot([0, 0], [0, 0], [0, line_length], color=line_color, linewidth=2, zorder=1)
    ax.text(-0.25, -line_length, -0.25, "Z", fontsize=12, color=line_color, zorder=1)

    ax.plot([0, 0], [0, -line_length], [0, 0], color=line_color, linewidth=2, zorder=1)
    ax.text(-0.07, 0, line_length + 0.1, "Y", fontsize=12, color=line_color, zorder=1)

    # Set axes limits.
    ax.set_xlim([0, 2.5])
    ax.set_ylim([-2.5, 0])
    ax.set_zlim([0, 2.5])

    return ax


def get_axes_2d() -> plt.Axes:
    """
    Get 2D axes.

    :return: 2D axes.
    :rtype: plt.Axes
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Remove axes.
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])

    # Plot 2D plot with two axes x,y as black arrows.
    ax.plot([0, 2.5], [0, 0], color="black", linewidth=2, zorder=1)
    ax.text(2.6, -0.1, "X", fontsize=12, color="black", zorder=1)

    ax.plot([0, 0], [0, 2.5], color="black", linewidth=2, zorder=1)
    ax.text(-0.1, 2.6, "Y", fontsize=12, color="black", zorder=1)

    # Set axes limits.
    ax.set_xlim([-2.5, 2.5])
    ax.set_ylim([-2.5, 2.5])

    return ax


def plot_points_3d(
    ax: Axes3D, points: ty.List[Point3D], color: ty.Union[ty.List[str], str]
) -> None:
    """
    Plot points in 3D.

    :param Axes3D ax: Axes to plot on.
    :param ty.List[Point3D] points: Points to plot.
    :param ty.Union[ty.List[str], str] color: Color of points.
    """
    points = np.array([[point.x, point.y, point.z] for point in points])
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color=color, s=0.01)


def plot_points_2d(
    ax: plt.Axes, points: ty.List[Point2D], color: ty.Union[ty.List[str], str], alpha: float = 1.0
) -> None:
    """
    Plot points in 2D.

    :param plt.Axes ax: Axes to plot on.
    :param ty.List[Point2D] points: Points to plot.
    :param ty.Union[ty.List[str], str] color: Color of points.
    """
    points = np.array([[point.x, point.y] for point in points])
    ax.scatter(points[:, 0], points[:, 1], color=color, s=0.01, alpha=alpha)


def save_plot_3d(ax: Axes3D, path: str, dpi: int, transparent_background: bool) -> None:
    """
    Save plot.

    :param Axes3D ax: Axes to plot on.
    :param str path: Path to save plot to.
    :param int dpi: DPI of output image.
    :param bool transparent_background: Transparent background.
    """
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=transparent_background)
    plt.clf()


def save_plot_2d(path: str, dpi: int, transparent_background: bool) -> None:
    """
    Save plot.

    :param str path: Path to save plot to.
    :param int dpi: DPI of output image.
    :param bool transparent_background: Transparent background.
    """
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=transparent_background)
    plt.clf()


# ==============================================================================
# Redefining these functions here because orientation axes system of matplotlib
# is different than for SVG.
# ==============================================================================


def get_points_on_surface_sphere(
    sphere: Sphere, num_phi: int, num_theta: int, filter_for_pov: bool = True
) -> ty.List[Point3D]:
    """
    Generate points on the surface of a sphere.

    :param Sphere sphere: The sphere.
    :param int num_phi: The resolution of the sphere in the phi direction.
    :param int num_theta: The resolution of the sphere in the theta direction.
    :param ty.Optional[Vector3D] filter_for_pov: If True, only return points
        that are on the surface of the sphere we can see from POV positive z-axis towards origin.
    :return: The points on the surface of the sphere.
    :rtype: ty.List[Point3D]
    """
    phis = [2 * math.pi * i / num_phi for i in range(num_phi + 1)]
    thetas = [math.pi * i / num_theta for i in range(num_theta + 1)]

    cx, cy, cz = sphere.center.x, sphere.center.y, sphere.center.z
    r = sphere.radius

    points = []
    for theta in thetas:
        for phi in phis:
            x = cx + r * math.sin(theta) * math.cos(phi)
            y = cy + r * math.sin(theta) * math.sin(phi)
            z = cz + r * math.cos(theta)

            # Only add points that are on the surface of the sphere we can see.
            if not filter_for_pov:
                points.append(Point3D(x, y, z))
                continue

            # Check if point is on the surface of the sphere we can see.
            if y <= sphere.center.y:
                points.append(Point3D(x, y, z))

    return points


def get_points_on_surface_cap(
    cap_type: CylinderCapType,
    center_cap: Point3D,
    radius_cap: float,
    normal_cap: Point3D,
    center_cylinder: Point3D,
    resolution: int,
) -> ty.List[Point3D]:
    """
    Generate points on the surface of the cap.

    :param CylinderCapType cap_type: The type of the cap.
    :param Point3D center_cap: The center of the cap.
    :param float radius_cap: The radius of the cap.
    :param Vector3D normal_cap: The normal vector of the cap.
    :param Point3D center_cylinder: The center of the cylinder.
    :param int resolution: The resolution of the cap.
    :return: The points on the surface of the cap.
    :rtype: ty.List[Point3D]
    """
    if cap_type == CylinderCapType.NO_CAP:
        return []

    elif cap_type == CylinderCapType.FLAT:
        circle = Circle3D(center_cap, radius_cap, normal_cap)
        return get_points_on_circumference_circle_3d(circle, resolution)

    elif cap_type == CylinderCapType.ROUND:
        sphere = Sphere(center_cap, radius_cap)
        plane = Plane3D(center_cap, normal_cap)
        points = get_points_on_surface_sphere(sphere, resolution, resolution, filter_for_pov=False)
        points = [
            point for point in points if not same_side_of_plane(plane, center_cylinder, point)
        ]
        return points

    else:
        raise ValueError(f"Unknown cap type: '{cap_type}'")


def get_points_on_surface_cylinder(cylinder: Cylinder, resolution: int) -> ty.List[Point3D]:
    """
    Generate points on the surface of the cylinder.

    :param Cylinder cylinder: The cylinder.
    :param int resolution: The resolution of the cylinder.
    :return: The points on the surface of the cylinder.
    :rtype: ty.List[Point3D]
    """
    normal = cylinder.end.create_vector(cylinder.start).normalize()
    centers = get_points_on_line_3d(Line3D(cylinder.start, cylinder.end), resolution)

    points = []
    for center in centers:
        circle = Circle3D(center, cylinder.radius, normal)
        points.extend(get_points_on_circumference_circle_3d(circle, resolution))

    # Get points on the caps.
    cap_resolution = max(int(resolution // 4), 2)
    cap_type = cylinder.cap_type
    cap_points = get_points_on_surface_cap(
        cap_type, cylinder.start, cylinder.radius, normal, cylinder.end, cap_resolution
    )
    points.extend(cap_points)
    cap_points = get_points_on_surface_cap(
        cap_type, cylinder.end, cylinder.radius, normal, cylinder.start, cap_resolution
    )
    points.extend(cap_points)

    return points


# ==============================================================================
# Steps of the algorithm explained visually
# ==============================================================================


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    pov_z = 15.0
    resolution = 100

    # Step 1.
    ax = get_axes_3d()

    cylinder1 = Cylinder(
        Point3D(0.0, -1.5, 0.3), Point3D(1.0, 0.0, 1.0), 0.2, CylinderCapType.ROUND
    )
    cylinder2 = Cylinder(
        Point3D(-0.1, -1.5, 0.3), Point3D(-1.0, 0.0, 1.5), 0.15, CylinderCapType.ROUND
    )
    sphere1 = Sphere(Point3D(1.0, 0.0, 1.0), 0.5)
    sphere2 = Sphere(Point3D(1.4, -0.2, 0.8), 0.3)

    points1 = get_points_on_surface_cylinder(cylinder1, resolution)
    points2 = get_points_on_surface_cylinder(cylinder2, resolution)
    points3 = get_points_on_surface_sphere(sphere1, resolution, resolution, filter_for_pov=False)
    points4 = get_points_on_surface_sphere(
        sphere2, resolution, int(resolution // 2.0), filter_for_pov=False
    )

    plot_points_3d(ax, points1, "black")
    plot_points_3d(ax, points2, "black")
    plot_points_3d(ax, points3, "black")
    plot_points_3d(ax, points4, "black")

    save_plot_3d(ax, os.path.join(args.o, "step1.png"), args.d, args.t)
    plt.clf()

    # Step 2.
    ax = get_axes_3d()

    filter_for_pov = True
    points1 = get_points_on_surface_cylinder(cylinder1, resolution)
    points2 = get_points_on_surface_cylinder(cylinder2, resolution)
    points3 = get_points_on_surface_sphere(
        sphere1, resolution, resolution, filter_for_pov=filter_for_pov
    )
    points4 = get_points_on_surface_sphere(
        sphere2, resolution, int(resolution // 2.0), filter_for_pov=filter_for_pov
    )

    points1 = [point for point in points1 if not point_is_inside_sphere(sphere1, point)]
    points2 = [point for point in points2 if not point_is_inside_cylinder(cylinder1, point)]
    points4 = [point for point in points4 if not point_is_inside_sphere(sphere1, point)]

    plot_points_3d(ax, points1, "black")
    plot_points_3d(ax, points2, "black")
    plot_points_3d(ax, points3, "black")
    plot_points_3d(ax, points4, "black")

    save_plot_3d(ax, os.path.join(args.o, "step2.png"), args.d, args.t)
    plt.clf()

    # Step 3.
    ax = get_axes_3d()

    def apply_perspective(points: ty.List[Point3D]) -> ty.List[Point3D]:
        new_points = []
        for point in points:
            s = pov_z / (pov_z - point.z)
            x = point.x * s
            y = point.y * s
            z = point.z * s
            new_points.append(Point3D(x, y, z))
        return new_points

    points1 = apply_perspective(points1)
    points2 = apply_perspective(points2)
    points3 = apply_perspective(points3)
    points4 = apply_perspective(points4)

    plot_points_3d(ax, points1, "black")
    plot_points_3d(ax, points2, "black")
    plot_points_3d(ax, points3, "black")
    plot_points_3d(ax, points4, "black")

    save_plot_3d(ax, os.path.join(args.o, "step3.png"), args.d, args.t)
    plt.clf()

    # Step 4.
    ax = get_axes_2d()

    plot_points_2d(ax, [Point2D(point.x, point.z) for point in points1], "black")
    plot_points_2d(ax, [Point2D(point.x, point.z) for point in points2], "black")
    plot_points_2d(ax, [Point2D(point.x, point.z) for point in points3], "black")
    plot_points_2d(ax, [Point2D(point.x, point.z) for point in points4], "black")

    save_plot_2d(os.path.join(args.o, "step4.png"), args.d, args.t)
    plt.clf()

    # Step 5.
    ax = get_axes_2d()

    points1 = [Point2D(point.x, point.z) for point in points1]
    points2 = [Point2D(point.x, point.z) for point in points2]
    points3 = [Point2D(point.x, point.z) for point in points3]
    points4 = [Point2D(point.x, point.z) for point in points4]

    # Calculate convex hull.
    hull1 = calculate_convex_hull(points1)
    points1 = [points1[i] for i in hull1]
    hull2 = calculate_convex_hull(points2)
    points2 = [points2[i] for i in hull2]
    hull3 = calculate_convex_hull(points3)
    points3 = [points3[i] for i in hull3]
    hull4 = calculate_convex_hull(points4)
    points4 = [points4[i] for i in hull4]

    # Draw hulls in order from furthest away to nearest.
    def draw_outline_hull(hull: ty.List[Point2D], color) -> None:
        for i in range(len(hull) - 1):
            ax.plot(
                [hull[i].x, hull[i + 1].x], [hull[i].y, hull[i + 1].y], color=color, linewidth=2
            )
        ax.plot([hull[-1].x, hull[0].x], [hull[-1].y, hull[0].y], color=color, linewidth=2)

    plot_points_2d(ax, points1, "black", alpha=0.2)
    plot_points_2d(ax, points2, "black", alpha=0.2)
    plot_points_2d(ax, points3, "black", alpha=0.2)
    plot_points_2d(ax, points4, "black", alpha=0.2)

    draw_outline_hull(points3, "black")
    draw_outline_hull(points4, "black")
    draw_outline_hull(points1, "black")
    draw_outline_hull(points2, "black")

    save_plot_2d(os.path.join(args.o, "step5.png"), args.d, args.t)
    plt.clf()

    # Step 6.
    ax = get_axes_2d()

    def draw_filled_hull(hull: ty.List[Point2D], color) -> None:
        x = [point.x for point in hull]
        y = [point.y for point in hull]
        ax.fill(x, y, color=color, alpha=1.0)

    draw_filled_hull(points3, "blue")
    draw_filled_hull(points4, "black")
    draw_filled_hull(points1, "red")
    draw_filled_hull(points2, "green")

    save_plot_2d(os.path.join(args.o, "step6"), args.d, args.t)


if __name__ == "__main__":
    main()
