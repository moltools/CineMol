# -*- coding: utf-8 -*-

"""Implements the quickhull algorithm for computing the convex hull of a set of points in 2D space.

Used as source: https://stackoverflow.com/questions/74812556/computing-quick-convex-hull-using-numba

Read more about the quickhull algorithm: https://surface.syr.edu/eecs_techreports/65/
"""

import typing as ty

from cinemol.geometry import Point2D


def argmin(vals: ty.List[float]) -> int:
    """Return the index of the minimum value in a list of floats.

    :param vals: The list of floats.
    :type vals: list[float]
    :return: The index of the minimum value in the list.
    :rtype: int
    """
    return min(range(len(vals)), key=lambda i: vals[i])


def argmax(vals: ty.List[float]) -> int:
    """Return the index of the maximum value in a list of floats.

    :param vals: The list of floats.
    :type vals: list[float]
    :return: The index of the maximum value in the list.
    :rtype: int
    """
    return max(range(len(vals)), key=lambda i: vals[i])


def arange(n: int) -> ty.List[int]:
    """Return a list of integers from 0 to n-1.

    :param n: The number of integers to return.
    :type n: int
    :return: The list of integers.
    :rtype: list[int]
    """
    return list(range(n))


def process(
    points: ty.List[Point2D],
    indices_of_points_to_consider: ty.List[int],
    index_a: int,
    index_b: int,
) -> ty.List[int]:
    """Recursively computes the convex hull of a set of points in 2D space.

    :param points: The list of points.
    :type points: list[Point2D]
    :param indices_of_points_to_consider: The indices of the points to consider.
    :type indices_of_points_to_consider: list[int]
    :param index_a: The index of the first point.
    :type index_a: int
    :param index_b: The index of the second point.
    :type index_b: int
    :return: The convex hull of the set of points.
    :rtype: list[int]
    """
    # Calculate the signed distance from each point in indices_of_points_to_consider
    # to the line between a and b.
    signed_dist = []
    for i in indices_of_points_to_consider:
        signed_dist.append(
            points[i]
            .subtract_point(points[index_a])
            .cross(points[index_b].subtract_point(points[index_a]))
        )

    # Find the points in indices_of_points_to_consider that are on the positive
    # side of the line between a and b.
    indices_on_positive_side = [
        i
        for s, i in zip(signed_dist, indices_of_points_to_consider)
        if (s > 0 and i != index_a and i != index_b)
    ]

    # If there are no points on the positive side of the line, return the line.
    if len(indices_on_positive_side) == 0:
        return [index_a, index_b]

    # Find the point in indices_of_points_to_consider that is farthest from the
    # line between a and b.
    index_c = max(zip(signed_dist, indices_of_points_to_consider))[1]

    # Recursively compute the convex hull of the points on the positive side of the line.
    return process(points, indices_on_positive_side, index_a, index_c)[:-1] + process(
        points, indices_on_positive_side, index_c, index_b
    )


def calculate_convex_hull(points: ty.List[Point2D]) -> ty.List[int]:
    """Calculate the convex hull of a set of points in 2D space.

    :param points: The list of points.
    :type points: list[Point2D]
    :return: The convex hull of the set of points.
    :rtype: list[int]
    """
    # Find the points with the minimum and maximum x-coordinates.
    min_index = argmin([p.x for p in points])
    max_index = argmax([p.x for p in points])

    # Recursively compute the convex hull of the points on the negative side of the line.
    return (
        process(points, arange(len(points)), min_index, max_index)[:-1]
        + process(points, arange(len(points)), max_index, min_index)[:-1]
    )
