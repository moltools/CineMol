#!/usr/bin/env python3
"""
Description:    Draw step-wise explanation of the algorithm.
Usage:          python3 draw_explanation_algorithm.py -h 
"""
import argparse 
import os
import typing as ty

import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 

from cinemol.geometry import Point2D, Point3D, Sphere, get_points_on_surface_sphere 
from cinemol.fitting import calculate_convex_hull

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
    parser.add_argument("-out", type=str, required=True, help="Path to output dir.")
    parser.add_argument("-dpi", type=int, required=False, default=300, help="DPI of output images.")
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
    ax = fig.add_subplot(111, projection='3d')

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
    ax.text(-0.1, 2.6, "Y", fontsize=12, color='black', zorder=1)

    # Set axes limits.
    ax.set_xlim([-2.5, 2.5])
    ax.set_ylim([-2.5, 2.5])

    return ax

def plot_sphere_surface_3d(ax: Axes3D, sphere: Sphere, color: str, alpha: float) -> None:
    """
    Plot sphere surface in 3D.

    :param Axes3D ax: Axes to plot on.
    :param Sphere sphere: Sphere to plot.
    :param str color: Color of sphere.
    :param float alpha: Alpha of sphere.
    """
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    cx, cy, cz = sphere.center.x, sphere.center.y, sphere.center.z
    radius = sphere.radius
    x = cx + radius * np.outer(np.cos(u), np.sin(v))
    y = cy + radius * np.outer(np.sin(u), np.sin(v))
    z = cz + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha)

def plot_points_3d(ax: Axes3D, points: ty.List[Point3D], color: ty.Union[ty.List[str], str]) -> None:
    """
    Plot points in 3D.
    
    :param Axes3D ax: Axes to plot on.
    :param ty.List[Point3D] points: Points to plot.
    :param ty.Union[ty.List[str], str] color: Color of points.
    """
    points = np.array([[point.x, point.y, point.z] for point in points])
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color=color, s=1)

def plot_points_2d(ax: plt.Axes, points: ty.List[Point3D], color: ty.Union[ty.List[str], str]) -> None:
    """
    Plot points in 2D.
    
    :param plt.Axes ax: Axes to plot on.
    :param ty.List[Point3D] points: Points to plot.
    :param ty.Union[ty.List[str], str] color: Color of points.
    """
    points = np.array([[point.x, point.z] for point in points])
    ax.scatter(points[:, 0], points[:, 1], color=color, s=1)

def save_plot_3d(ax: Axes3D, path: str, dpi: int) -> None:
    """
    Save plot.
    
    :param Axes3D ax: Axes to plot on.
    :param str path: Path to save plot to.
    :param int dpi: DPI of output image.
    """
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

def save_plot_2d(path: str, dpi: int) -> None:
    """
    Save plot.
    
    :param str path: Path to save plot to.
    :param int dpi: DPI of output image.
    """
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

# ==============================================================================
# Steps of the algorithm explained visually
# ==============================================================================

def step1(path: str, dpi: int) -> None:
    """
    Step 1 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    """
    ax = get_axes_3d()
    plot_sphere_surface_3d(ax, Sphere(Point3D(0.0, 0.0, 0.0), 0.75), "blue", 0.25)
    plot_sphere_surface_3d(ax, Sphere(Point3D(0.25, 0.0, 0.5), 0.75), "red", 0.25)
    save_plot_3d(ax, path, dpi)

def step2(path: str, dpi: int) -> ty.List[Point3D]:
    """
    Step 2 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :return: Visible points.
    :rtype: ty.List[Point3D]
    """
    ax = get_axes_3d()
    sphere1 = Sphere(Point3D(0.0, 0.0, 0.0), 0.75)
    sphere2 = Sphere(Point3D(0.25, 0.0, 0.5), 0.75)
    plot_sphere_surface_3d(ax, sphere1, "blue", 0.25)
    plot_sphere_surface_3d(ax, sphere2, "red", 0.25)

    num_phis, num_thetas = 15, 15
    points = get_points_on_surface_sphere(sphere1, num_phis, num_thetas, filter_for_pov=False)
    points = [p for p in points if p.y <= sphere1.center.y]
    visible_points, colors = [], []
    for point in points:
        dist = point.calculate_distance(sphere2.center)
        if dist <= sphere2.radius:
            colors.append("red")
        else:
            colors.append("black")
            visible_points.append(point)

    plot_points_3d(ax, points, colors)
    save_plot_3d(ax, path, dpi)
    return visible_points

def step3(path: str, dpi: int, visible_points: ty.List[Point3D]) -> None:
    """
    Step 3 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :param ty.List[Point3D] visible_points: Visible points.
    """
    ax = get_axes_2d()
    plot_points_2d(ax, visible_points, "black")
    save_plot_2d(path, dpi)

def step4(path: str, dpi: int, visible_points: ty.List[Point3D]) -> ty.List[Point2D]:
    """
    Step 4 of the algorithm.

    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :param ty.List[Point3D] visible_points: Visible points.
    :return: Convex hull.
    :rtype: ty.List[Point2D]
    """
    ax = get_axes_2d()
    plot_points_2d(ax, visible_points, "black")

    # Calculate convex hull.
    points = [Point2D(point.x, point.z) for point in visible_points]
    convex_hull = calculate_convex_hull(points)
    hull = [points[i] for i in convex_hull]

    # Add hull as line segments.
    for i in range(len(hull) - 1):
        ax.plot([hull[i].x, hull[i + 1].x], [hull[i].y, hull[i + 1].y], color="blue", linewidth=2)
    ax.plot([hull[-1].x, hull[0].x], [hull[-1].y, hull[0].y], color="blue", linewidth=2)

    save_plot_2d(path, dpi) 

    return hull 

def step5(path: str, dpi: int, hull: ty.List[Point2D]) -> None:
    """
    Step 5 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :param ty.List[Point2D] hull: Hull.
    """
    ax = get_axes_2d()

    # Plot sphere 2 as circle on x/y axis.
    sphere = Sphere(Point3D(0.25, 0.0, 0.5), 0.75)
    circle = plt.Circle((sphere.center.x, sphere.center.z), sphere.radius, color="red", fill=True)
    ax.add_artist(circle)

    # Plot hull as polygon on top of sphere.
    x = [point.x for point in hull]
    y = [point.y for point in hull]
    ax.fill(x, y, color="blue", alpha=1.0)

    save_plot_2d(path, dpi)

def main() -> None:
    """
    Main function.
    """
    args = cli()
    os.makedirs(args.out, exist_ok=True)

    # Step 1.
    step1(os.path.join(args.out, "step1.png"), args.dpi)

    # Step 2.
    visible_points = step2(os.path.join(args.out, "step2.png"), args.dpi)

    # Step 3.
    step3(os.path.join(args.out, "step3.png"), args.dpi, visible_points)

    # Step 4.
    hull = step4(os.path.join(args.out, "step4.png"), args.dpi, visible_points)

    # Step 5.
    step5(os.path.join(args.out, "step5.png"), args.dpi, hull)

    exit(0)

if __name__ == "__main__":
    main()