import typing as ty 

from cinemol.geometry import Point2D 

def process(S: ty.List[Point2D], P: ty.List[int], a: int, b: int) -> ty.List[int]:
    signed_dist = []
    for i in P: 
        p = S[i]
        p = Point2D(p.x - S[a].x, p.y - S[a].y) # Subtract a from p 
        p = Point2D(p.x * S[b].y, p.y * S[b].x) # Cross product with b
        p = Point2D(p.x - S[a].x, p.y - S[a].y) # Subtract a from p
        signed_dist.append(p)

    K = [i for s, i in zip(signed_dist, P) if s > 0 and i != a and i != b]

    if len(K) == 0: 
        return (a, b)
    
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