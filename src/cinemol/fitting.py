"""
Quickhull algorithm for computing the convex hull of a set of points in 2D space.

Source: https://stackoverflow.com/questions/74812556/computing-quick-convex-hull-using-numba
"""
import typing as ty 

from cinemol.geometry import Point2D  

def argmin(vals: ty.List[float]) -> int:
    """
    Returns the index of the minimum value in a list of floats.
    
    :param list[float] vals: The list of floats.
    :return: The index of the minimum value in the list.
    :rtype: int
    """
    return min(range(len(vals)), key=lambda i: vals[i])

def argmax(vals: ty.List[float]) -> int:
    """
    Returns the index of the maximum value in a list of floats.
    
    :param list[float] vals: The list of floats.
    :return: The index of the maximum value in the list.
    :rtype: int
    """
    return max(range(len(vals)), key=lambda i: vals[i])

def arange(n: int) -> ty.List[int]:
    """
    Returns a list of integers from 0 to n-1.
    
    :param int n: The number of integers to return.
    :return: A list of integers from 0 to n-1.
    :rtype: list[int]
    """
    return list(range(n))

def process(S: ty.List[Point2D], P: ty.List[int], a: int, b: int) -> ty.List[int]:
    """
    Recursively computes the convex hull of a set of points in 2D space.
    
    :param list[Point2D] S: The list of points.
    :param list[int] P: The list of indices of points to consider.
    :param int a: The index of the first point in the set.
    :param int b: The index of the second point in the set.
    :return: The convex hull of the set of points.
    :rtype: list[int]
    """
    # Calculate the signed distance from each point in P to the line between a and b.
    signed_dist = []
    for i in P: 
        signed_dist.append(S[i].subtract_point(S[a]).cross(S[b].subtract_point(S[a])))

    # Find the points in P that are on the positive side of the line between a and b.
    K = [i for s, i in zip(signed_dist, P) if s > 0 and i != a and i != b]

    # If there are no points on the positive side of the line, return the line.
    if len(K) == 0: 
        return (a, b)
    
    # Find the point in P that is farthest from the line between a and b.
    c = max(zip(signed_dist, P))[1]

    # Recursively compute the convex hull of the points on the positive side of the line.
    return process(S, K, a, c)[:-1] + process(S, K, c, b)

def calculate_convex_hull(S) -> ty.List[int]:
    """
    Calculates the convex hull of a set of points in 2D space.
    
    :param list[Point2D] S: The list of points.
    :return: The convex hull of the set of points.
    :rtype: list[int]
    """
    # Find the points with the minimum and maximum x-coordinates.
    a = argmin([p.x for p in S])

    # Find the point with the maximum x-coordinate. 
    max_index = argmax([p.x for p in S])

    # Recursively compute the convex hull of the points on the negative side of the line.
    return (
        process(S, arange(len(S)), a, max_index)[:-1] + 
        process(S, arange(len(S)), max_index, a)[:-1]
    )