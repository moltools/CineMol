"""
Source: https://stackoverflow.com/questions/74812556/computing-quick-convex-hull-using-numba
"""
import typing as ty 

import numpy as np

def process(S: np.ndarray, P: np.ndarray, a: int, b: int) -> ty.Tuple[int, ...]:
    signed_dist = np.cross(S[P] - S[a], S[b] - S[a])
    K = [i for s, i in zip(signed_dist, P) if s > 0 and i != a and i != b]

    if len(K) == 0:
        return (a, b)

    c = max(zip(signed_dist, P))[1]
    return process(S, K, a, c)[:-1] + process(S, K, c, b)

def quick_hull_2d(S: np.ndarray) -> np.ndarray:
    a = np.argmin(S[:,0])
    max_index = np.argmax(S[:,0])

    return (
        process(S, np.arange(S.shape[0]), a, max_index)[:-1] + 
        process(S, np.arange(S.shape[0]), max_index, a)[:-1]
    )