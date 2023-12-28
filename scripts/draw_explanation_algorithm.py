#!/usr/bin/env python3
import argparse 
import typing as ty

import matplotlib.pyplot as plt 
import numpy as np 

def cli() -> argparse.Namespace:
    """
    Command line interface for this script.
    
    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-out", type=str, required=True, help="Path to output dir.")
    parser.add_argument("-dpi", type=int, default=300, help="DPI of output images.")
    return parser.parse_args()

def step1(path: str, dpi: int) -> None:
    """
    Step 1 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    line_length = 2.0
    line_color = "black"
    ax.plot([0, line_length], [0, 0], [0, 0], color=line_color, linewidth=2)
    ax.text(line_length + 0.1, -0.1, 0, "X", fontsize=12, color=line_color)
    ax.plot([0, 0], [0, 0], [0, line_length], color=line_color, linewidth=2)
    ax.text(-0.25, -line_length, -0.25, "Z", fontsize=12, color=line_color)
    ax.plot([0, 0], [0, -line_length], [0, 0], color=line_color, linewidth=2)
    ax.text(-0.07, 0, line_length + 0.1, "Y", fontsize=12, color=line_color)
    ax.set_xlim([0, 2.5])
    ax.set_ylim([-2.5, 0])
    ax.set_zlim([0, 2.5])

    # Add intersecting spheres (blue and red)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Sphere 1
    radius = 0.75
    x_center, y_center, z_center = 0, 0, 0
    x = x_center + radius * np.outer(np.cos(u), np.sin(v))
    y = y_center + radius * np.outer(np.sin(u), np.sin(v))
    z = z_center + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color="blue", alpha=0.25)

    # Sphere 2
    radius = 0.75
    x_center, y_center, z_center = 0.25, 0, 0.5
    x = x_center + radius * np.outer(np.cos(u), np.sin(v))
    y = y_center + radius * np.outer(np.sin(u), np.sin(v))
    z = z_center + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color="red", alpha=0.25)

    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

def step2(path: str, dpi: int) -> np.ndarray:
    """
    Step 2 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :return: Visible points.
    :rtype: np.ndarray
    """
    # Plot 3D plot with three axes x,y,z as black arrows.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    line_length = 2.0
    line_color = "black"
    ax.plot([0, line_length], [0, 0], [0, 0], color=line_color, linewidth=2)
    ax.text(line_length + 0.1, -0.1, 0, "X", fontsize=12, color=line_color)
    ax.plot([0, 0], [0, 0], [0, line_length], color=line_color, linewidth=2)
    ax.text(-0.25, -line_length, -0.25, "Z", fontsize=12, color=line_color)
    ax.plot([0, 0], [0, -line_length], [0, 0], color=line_color, linewidth=2)
    ax.text(-0.07, 0, line_length + 0.1, "Y", fontsize=12, color=line_color)
    ax.set_xlim([0, 2.5])
    ax.set_ylim([-2.5, 0])
    ax.set_zlim([0, 2.5])

    # Add intersecting spheres (blue and red)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Sphere 1
    sphere1_radius = 0.75
    sphere1_x_center, sphere1_y_center, sphere1_z_center = 0, 0, 0
    x = sphere1_x_center + sphere1_radius * np.outer(np.cos(u), np.sin(v))
    y = sphere1_y_center + sphere1_radius * np.outer(np.sin(u), np.sin(v))
    z = sphere1_z_center + sphere1_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color="blue", alpha=0.25)

    # Sphere 2
    sphere2_radius = 0.75
    sphere2_x_center, sphere2_y_center, sphere2_z_center = 0.25, 0, 0.5
    x = sphere2_x_center + sphere2_radius * np.outer(np.cos(u), np.sin(v))
    y = sphere2_y_center + sphere2_radius * np.outer(np.sin(u), np.sin(v))
    z = sphere2_z_center + sphere2_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color="red", alpha=0.25)

    # Generate points on the surface of the blue sphere based on polar and azimuthal angles
    num_points = 20
    phi = np.linspace(0, 2 * np.pi, num_points)
    theta = np.linspace(0, np.pi, num_points)
    phi, theta = np.meshgrid(phi, theta)
    x_points = sphere1_x_center + sphere1_radius * np.sin(theta) * np.cos(phi)
    y_points = sphere1_y_center + sphere1_radius * np.sin(theta) * np.sin(phi)
    z_points = sphere1_z_center + sphere1_radius * np.cos(theta)
    points = np.array([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T

    # Only keep points that are visible to viewer.
    points = points[points[:, 1] <= sphere1_z_center]

    visible_points = []
    colors = []
    for point in points:
        dist_to_sphere2 = np.sqrt((point[0] - sphere2_x_center)**2 + (point[1] - sphere2_y_center)**2 + (point[2] - sphere2_z_center)**2)
        if dist_to_sphere2 <= sphere2_radius:
            colors.append("red")
        else:
            colors.append("black")
            visible_points.append(point)

    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color=colors, s=1)

    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

    visible_points = np.array(visible_points)

    return visible_points

def step3(path: str, dpi: int, visible_points: np.ndarray) -> None:
    """
    Step 3 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :param np.ndarray visible_points: Visible points.
    """
    # Plot visible points on 2D.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-2.5, 2.5])
    ax.set_ylim([-2.5, 2.5])

    # Line along x axis
    ax.plot([0, 2.5], [0, 0], color="black", linewidth=2)
    ax.text(2.6, -0.1, "X", fontsize=12, color="black")

    # Line along y axis
    ax.plot([0, 0], [0, 2.5], color="black", linewidth=2)
    ax.text(-0.1, 2.6, "Y", fontsize=12, color='black')

    ax.scatter(visible_points[:, 0], visible_points[:, 2], color="black", s=1)
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

def step4(path: str, dpi: int, visible_points: np.ndarray) -> None:
    """
    Step 4 of the algorithm.

    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    :param np.ndarray visible_points: Visible points.
    """
    # Plot visible points on 2D.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-2.5, 2.5])
    ax.set_ylim([-2.5, 2.5])

    # Line along x axis
    ax.plot([0, 2.5], [0, 0], color="black", linewidth=2)
    ax.text(2.6, -0.1, "X", fontsize=12, color="black")

    # Line along y axis
    ax.plot([0, 0], [0, 2.5], color="black", linewidth=2)
    ax.text(-0.1, 2.6, "Y", fontsize=12, color="black")

    ax.scatter(visible_points[:, 0], visible_points[:, 2], color="black", s=1)

    class Point2D:
        def __init__(self, x: float, y: float) -> None:
            self.x = x 
            self.y = y

    def process(S: ty.List[Point2D], P: ty.List[int], a: int, b: int) -> ty.List[int]:
        signed_dist = []
        for i in P: signed_dist.append(S[i].subtract_point(S[a]).cross(S[b].subtract_point(S[a])))
        K = [i for s, i in zip(signed_dist, P) if s > 0 and i != a and i != b]
        if len(K) == 0: return (a, b)
        c = max(zip(signed_dist, P))[1]
        return process(S, K, a, c)[:-1] + process(S, K, c, b)

    def argmin(vals: ty.List[float]) -> int:
        return min(range(len(vals)), key=lambda i: vals[i])

    def argmax(vals: ty.List[float]) -> int:
        return max(range(len(vals)), key=lambda i: vals[i])

    def arange(n: int) -> ty.List[int]:
        return list(range(n))

    def quickhull_2d(S) -> ty.List[int]:
        a = argmin([p.x for p in S])
        max_index = argmax([p.x for p in S])
        return (
            process(S, arange(len(S)), a, max_index)[:-1] + 
            process(S, arange(len(S)), max_index, a)[:-1]
        )

    convex_hull = quickhull_2d([Point2D(point[0], point[2]) for point in visible_points])
    hull = [visible_points[i] for i in convex_hull]

    # Add hull as line segments.
    for i in range(len(hull) - 1):
        ax.plot([hull[i][0], hull[i + 1][0]], [hull[i][2], hull[i + 1][2]], color="blue", linewidth=2)
    ax.plot([hull[-1][0], hull[0][0]], [hull[-1][2], hull[0][2]], color="blue", linewidth=2)

    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

def step5(path: str, dpi: int) -> None:
    """
    Step 5 of the algorithm.
    
    :param str path: Path to save image to.
    :param int dpi: DPI of output image.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-2.5, 2.5])
    ax.set_ylim([-2.5, 2.5])

    # Line along x axis
    ax.plot([0.83, 2.5], [0, 0], color="black", linewidth=2)
    ax.text(2.6, -0.1, "X", fontsize=12, color="black")

    # Line along y axis
    ax.plot([0, 0], [1.2, 2.5], color='black', linewidth=2)
    ax.text(-0.1, 2.6, "Y", fontsize=12, color="black")

    # Plot sphere 2 as circle on x/y axis.
    sphere2_radius = 0.75
    sphere2_x_center, sphere2_y_center, sphere2_z_center = 0.25, 0, 0.5
    circle = plt.Circle((sphere2_x_center, sphere2_z_center), sphere2_radius, color="red", fill=True)
    ax.add_artist(circle)

    # Plot hull as polygon.
    hull = np.array(hull)
    ax.fill(hull[:, 0], hull[:, 2], color="blue", alpha=1.0)

    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, transparent=True)
    plt.clf()

def main() -> None:
    args = cli()
    exit(0)

if __name__ == "__main__":
    main()