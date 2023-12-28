"""
Quickhull algorithm for computing the convex hull of a set of points in 2D space.

Source: https://stackoverflow.com/questions/74812556/computing-quick-convex-hull-using-numba
"""
import typing as ty 

import numpy as np

def process(S: np.ndarray, P: np.ndarray, a: int, b: int) -> ty.Tuple[int, ...]:
    """
    Process the points in P that are on the right side of the line from a to b.
    
    :param np.ndarray S: The set of points.
    :param np.ndarray P: The set of indices of points to process.
    :param int a: The index of the first point of the line.
    :param int b: The index of the second point of the line.
    :return: The indices of the points in P that are on the right side of the line from a to b.
    :rtype: ty.Tuple[int, ...]
    """
    # Find the point c that is furthest from the line from a to b.
    signed_dist = np.cross(S[P] - S[a], S[b] - S[a])
    K = [i for s, i in zip(signed_dist, P) if s > 0 and i != a and i != b]

    # If there are no points on the right side, return the line from a to b.
    if len(K) == 0:
        return (a, b)

    # Find the point c that is furthest from the line from a to b.
    c = max(zip(signed_dist, P))[1]

    # Process the points on the right side of the line from a to c.
    return process(S, K, a, c)[:-1] + process(S, K, c, b)

def quick_hull_2d(S: np.ndarray) -> np.ndarray:
    """
    Compute the convex hull of a set of points in 2D space.
    
    :param np.ndarray S: The set of points.
    :return: The indices of the points that form the convex hull.
    :rtype: np.ndarray
    """
    # Find the indices of the leftmost and rightmost points.
    a = np.argmin(S[:,0])
    max_index = np.argmax(S[:,0])

    # Process the points on the right side of the line from a to max_index.
    return (
        process(S, np.arange(S.shape[0]), a, max_index)[:-1] + 
        process(S, np.arange(S.shape[0]), max_index, a)[:-1]
    )