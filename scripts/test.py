import sys

import rdkit 
from rdkit import Chem

import typing as ty
from dataclasses import dataclass
from enum import Enum

import numpy as np

@dataclass 

class Color:
    r: int
    g: int
    b: int

    def __str__(self) -> str:
        return f"rgb({self.r}, {self.g}, {self.b})"

@dataclass
class Point3D:
    x: float
    y: float
    z: float

@dataclass 
class Vector3D:
    x: float
    y: float
    z: float

@dataclass
class Point2D:
    x: float
    y: float

    def cross(self, other: "Point2D") -> float:
        return self.x * other.y - self.y * other.x
    
    def subtract_point(self, other: "Point2D") -> "Point2D":
        return Point2D(self.x - other.x, self.y - other.y)

@dataclass 
class Node:
    index: int 
    center: Point3D
    radius: float = 0.5
    color: Color = Color(0, 0, 0)

    def intersects(self, other: "Node") -> bool:
        distance = ((self.center.x - other.center.x) ** 2 + (self.center.y - other.center.y) ** 2) ** 0.5
        return distance <= self.radius + other.radius
    
    def encapsulates(self, other: "Node") -> bool:
        distance = ((self.center.x - other.center.x) ** 2 + (self.center.y - other.center.y) ** 2) ** 0.5
        return distance + other.radius <= self.radius
    
@dataclass 
class Edge:
    index: int 
    start_index: int
    start: Point3D
    end_index: int 
    end: Point3D
    radius: float = 0.2
    color: Color = Color(0, 0, 0)

@dataclass
class ViewBox:
    min_x: float 
    min_y: float 
    width: float 
    height: float

    def __str__(self) -> str:
        return f"viewBox=\"{self.min_x:.3f} {self.min_y:.3f} {self.width:.3f} {self.height:.3f}\""

@dataclass 
class Svg:
    view_box: ViewBox  
    version: float = 1.0
    encoding: str = "UTF-8"

    def header(self) -> str:
        return (
            f"<?xml version=\"{self.version}\" encoding=\"{self.encoding}\"?>\n"
            f"<svg xmlns=\"http://www.w3.org/2000/svg\" {self.view_box}>"
        )
    
    def footer(self) -> str:
        return "</svg>"
    
def generate_points_on_sphere(center: Point3D, radius: float, resolution: int) -> ty.List[Point3D]:
    num_points = resolution
    sphere1_x_center = center.x
    sphere1_y_center = center.y
    sphere1_z_center = center.z
    phi = np.linspace(0, 2 * np.pi, num_points)
    theta = np.linspace(0, np.pi, num_points)
    phi, theta = np.meshgrid(phi, theta)
    x_points = sphere1_x_center + radius * np.sin(theta) * np.cos(phi)
    y_points = sphere1_y_center + radius * np.sin(theta) * np.sin(phi)
    z_points = sphere1_z_center + radius * np.cos(theta)
    points = np.array([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T
    points = [Point3D(point[0], point[1], point[2]) for point in points]
    return points

def point_inside_sphere(point: Point3D, sphere_center: Point3D, sphere_radius: float) -> bool:
    dist = ((point.x - sphere_center.x) ** 2 + (point.y - sphere_center.y) ** 2 + (point.z - sphere_center.z) ** 2) ** 0.5
    return dist <= sphere_radius

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
    # return process(S, arange(len(S)), a, max_index)[:-1] + process(S, arange(len(S)), max_index, a)[:-1]
    return (
        process(S, arange(len(S)), a, max_index)[:-1] + 
        process(S, arange(len(S)), max_index, a)[:-1]
    )

    
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

def create_polygon(node: Node, other_nodes: ty.List[Node]) -> ty.List[Point2D]:
    points = generate_points_on_sphere(node.center, node.radius, 100)

    if len(other_nodes) == 0:
        return points

    visible_points = []
    for point in points:
        if all([not point_inside_sphere(point, other_node.center, other_node.radius) for other_node in other_nodes]):
            visible_points.append(point)

    if len(visible_points) == 0:
        return []

    points_2d = [Point2D(point.x, point.y) for point in visible_points]
    points_2d = quickhull_2d(points_2d)

    return [visible_points[i] for i in points_2d]

class Cap(Enum):
    Round = 1
    Flat = 2
    NoCap = 3

@dataclass
class Circle3D:
    center: Point3D
    radius: float
    normal: Vector3D

def normalize(V):
    return V / np.linalg.norm(V)

def gram_schmidt(N):
    V = np.random.rand(3)
    V = V - np.dot(V, N) * N # Create V orthogonal to N
    V = normalize(V)
    W = np.cross(N, V) # Create W orthogonal to N and V
    return V, W

def calculate_normal(p1: Point3D, p2: Point3D) -> Vector3D:
    return Vector3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)

def calculate_points_along_line(start: Point3D, end: Point3D, num_points: int) -> ty.List[Point3D]:
    x = np.linspace(start.x, end.x, num_points)
    y = np.linspace(start.y, end.y, num_points)
    z = np.linspace(start.z, end.z, num_points)
    points = [Point3D(x[i], y[i], z[i]) for i in range(num_points)]
    return points

def generate_points_on_circle(circle: Circle3D, num_points: int):
    normal_vector = np.array([circle.normal.x, circle.normal.y, circle.normal.z])
    center = np.array([circle.center.x, circle.center.y, circle.center.z])
    radius = circle.radius
    N = normalize(normal_vector)
    V, W = gram_schmidt(N) # TODO: for cylinder we calculate this now many times...
    angles = np.linspace(0, 2 * np.pi, num_points)
    points = center + radius * np.outer(np.cos(angles), V) + radius * np.outer(np.sin(angles), W)
    points = [Point3D(point[0], point[1], point[2]) for point in points]
    return points

def point_side_of_plane(plane_normal, plane_point, point1, point2):
    vector1 = np.array(point1) - np.array(plane_point)
    vector2 = np.array(point2) - np.array(plane_point)
    dot_product1 = np.dot(vector1, plane_normal)
    dot_product2 = np.dot(vector2, plane_normal)
    return np.sign(dot_product1) == np.sign(dot_product2)

def calc_cap(normal: Vector3D, center: Point3D, other: Point3D, radius: float, resolution: int, cap: Cap, points: ty.List[Point3D]):
    normal = np.array([normal.x, normal.y, normal.z])
    match cap:
        case Cap.Round:
            points = generate_points_on_sphere(center, radius, resolution)
            points = np.array([[point.x, point.y, point.z] for point in points])
            plane_normal = np.array([other.x - center.x, other.y - center.y, other.z - center.z])
            plane_normal = normalize(plane_normal)
            center = np.array([center.x, center.y, center.z])
            other = np.array([other.x, other.y, other.z])
            points = [point for point in points if not point_side_of_plane(plane_normal, center, point, other)]
            # TODO: can probably speed up calc above by just calculating two opositie points on sphere and then just checking if point is on the same side as those two points
            points = [Point3D(point[0], point[1], point[2]) for point in points]
        case Cap.Flat:
            points = []
            radii = np.linspace(0, radius, resolution)
            for radius in radii:
                cap1 = Circle3D(center, radius, normal)
                points.extend(generate_points_on_circle(cap1, resolution))
        case Cap.NoCap:
            points= []
    return points 

def generate_points_on_cylinder(start: Point3D, end: Point3D, radius: float, resolution: int, start_cap: Cap, end_cap: Cap) -> ty.List[Point3D]:
    normal = calculate_normal(start, end)
    centers = calculate_points_along_line(start, end, resolution)
    points = []
    for center in centers:
        circle = Circle3D(center, radius, normal)
        points_on_circle = generate_points_on_circle(circle, resolution)
        points.extend(points_on_circle)
    
    start_cap_points = calc_cap(normal, start, end, radius, resolution, start_cap, points)
    end_cap_points = calc_cap(normal, end, start, radius, resolution, end_cap, points)
    points.extend(start_cap_points)
    points.extend(end_cap_points)
    return points

class Scene:
    def __init__(self, view_box: ViewBox = None):
        self.nodes = []
        self.edges = []

    def add_node(self, node: Node):
        self.nodes.append(node)

    def add_edge(self, edge: Edge):
        self.edges.append(edge)

    def circle_to_svg(self, center: Point2D, radius: float, color: Color) -> str:
        return f"<circle cx=\"{center.x:.3f}\" cy=\"{center.y:.3f}\" r=\"{radius:.3f}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.01\"/>"

    def polygon_to_svg(self, points: ty.List[Point2D], color: Color) -> str:
        points_str = " ".join([f"{point.x:.3f},{point.y:.3f}" for point in points])
        return f"<polygon points=\"{points_str}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.01\"/>"
        
    def _sort_nodes(self) -> None:
        return sorted(self.nodes, key=lambda node: node.center.z, reverse=False)

    def _calculate_view_box(self) -> ViewBox:
        x = [node.center.x for node in self.nodes]
        y = [node.center.y for node in self.nodes]
        min_x, min_y, max_x, max_y = min(x), min(y), max(x), max(y)
        margin = 5
        min_x -= margin
        min_y -= margin
        max_x += margin
        max_y += margin
        return ViewBox(min_x, min_y, max_x - min_x, max_y - min_y)
    
    def _generate_style(self) -> str:
        return "" # TODO
    
    def _generate_defs(self) -> str:
        return "" # TODO
    
    def process_bonds(self) -> ty.List[str]:
        # Filter out edges that are not connected to the current node, and edges that are connected to the previous nodes.
        # prev_node_inds = [prev_node.index for prev_node in prev_nodes]
        # edges = self.edges 
        # edges = [
        #     edge for edge 
        #     in edges 
        #     if (edge.start_index == node.index or edge.end_index == node.index) and
        #         edge.start_index not in prev_node_inds and edge.end_index not in prev_node_inds
            
        # ]

        edge_strs = []

        for i, edge in enumerate(self.edges):
            # start_node = [node for node in self.nodes if node.index == edge.start_index][0]
            # end_node = [node for node in self.nodes if node.index == edge.end_index][0]
            start = edge.start 
            start_atom = [node for node in self.nodes if node.index == edge.start_index][0]
            start_color = start_atom.color
            end = edge.end
            end_atom = [node for node in self.nodes if node.index == edge.end_index][0]
            end_color = end_atom.color

            # Get start/end based on which one is furthest away (most negative z value) 
            mid = Point3D((start.x + end.x) / 2, (start.y + end.y) / 2, (start.z + end.z) / 2)

            # first part of bond
            points = generate_points_on_cylinder(start, mid, 0.1, 50, Cap.NoCap, Cap.NoCap)
            filtered_points = []
            for point in points:
                if any([point_inside_sphere(point, node.center, node.radius) for node in self.nodes]):
                    continue
                filtered_points.append(point)
            hull = quickhull_2d([Point2D(point.x, point.y) for point in filtered_points])
            filtered_points = [filtered_points[i] for i in hull]
            if len(filtered_points) == 0: continue
            edge_str = self.polygon_to_svg(filtered_points, start_color)
            start_mid_mid = Point3D((start.x + mid.x) / 2, (start.y + mid.y) / 2, (start.z + mid.z) / 2)
            edge_strs.append((edge_str, start_mid_mid))

            # second part of bond
            points = generate_points_on_cylinder(mid, end, 0.1, 50, Cap.NoCap, Cap.NoCap)
            filtered_points = []
            for point in points:
                if any([point_inside_sphere(point, node.center, node.radius) for node in self.nodes]):
                    continue
                filtered_points.append(point)
            hull = quickhull_2d([Point2D(point.x, point.y) for point in filtered_points])
            filtered_points = [filtered_points[i] for i in hull]
            if len(filtered_points) == 0: continue
            edge_str = self.polygon_to_svg(filtered_points, end_color)
            mid_end_mid = Point3D((mid.x + end.x) / 2, (mid.y + end.y) / 2, (mid.z + end.z) / 2)
            edge_strs.append((edge_str, mid_end_mid))

            print(f"{i}".zfill(5), end="\r")

        return edge_strs # with sort pos
    
    def _generate_objs(self) -> str:
        nodes = self._sort_nodes()
        
        svgs = [] # TODO: should keep intermediate objects instead of already the svgs.
        for i, node in enumerate(nodes):

            if i == 0:
                svgs.append((self.circle_to_svg(Point2D(node.center.x, node.center.y), node.radius, node.color), node.center))
                continue
            
            prev_nodes = nodes[:i]
            if any([prev_node.encapsulates(node) for prev_node in prev_nodes]):
                continue

            prev_nodes = [prev_node for prev_node in prev_nodes if node.intersects(prev_node)]
            if len(prev_nodes) == 0:
                svgs.append((self.circle_to_svg(Point2D(node.center.x, node.center.y), node.radius, node.color), node.center))
                continue
                
            polygon = create_polygon(node, prev_nodes)
            if len(polygon) == 0:
                continue
                
            svgs.append((self.polygon_to_svg(polygon, node.color), node.center))

            print(f"{i}".zfill(5), end="\r")
        
        bond_svgs = self.process_bonds()
        svgs.extend(bond_svgs)

        # sort by z value
        svgs = sorted(svgs, key=lambda svg: svg[1].z, reverse=False)
        svgs = [svg[0] for svg in svgs]

        return "\n".join(svgs)

    def to_svg(self) -> str:
        view_box = self._calculate_view_box()
        svg = Svg(view_box=view_box)
        svg_str = "{header}\n<defs>\n<style>\n{style}\n</style>\n{defs}\n</defs>\n{objs}\n{footer}"
        svg_str = svg_str.format(
            header=svg.header(),
            style=self._generate_style(),
            defs=self._generate_defs(),
            objs=self._generate_objs(),
            footer=svg.footer()
        )
        return svg_str

file_path = sys.argv[1]
mol = Chem.MolFromMolFile(file_path, removeHs=False)
pos = mol.GetConformer().GetPositions()

def atom_symbol_to_radius(atom_symbol: str) -> float:
    return {
        "H": 0.6,
        "C": 1.20,
        "N": 1.20,
        "O": 0.73,
        "F": 0.71,
        "P": 1.06,
        "S": 1.40,
        "Cl": 0.99,
        "Br": 1.14,
        "I": 1.33,
    }[atom_symbol]

def atom_symbol_to_color(atom_symbol: str) -> Color:
    return {
        "H": Color(255, 255, 255),
        "C": Color(200, 200, 200),
        "N": Color(143, 143, 255),
        "O": Color(255, 0, 0),
        "F": Color(144, 255, 144),
        "P": Color(255, 165, 0),
        "S": Color(255, 255, 0),
        "Cl": Color(0, 255, 0),
        "Br": Color(165, 42, 42),
        "I": Color(148, 0, 211),
    }[atom_symbol]

scene = Scene()

for i, atom in enumerate(mol.GetAtoms()):
    x, y, z = pos[i]
    atom_symbol = atom.GetSymbol()
    radius = atom_symbol_to_radius(atom_symbol) / 3
    color = atom_symbol_to_color(atom_symbol)
    node = Node(atom.GetIdx(), Point3D(x, y, z), radius, color)
    scene.add_node(node)

for i, bond in enumerate(mol.GetBonds()):
    start_index = bond.GetBeginAtomIdx()
    start_pos = pos[bond.GetBeginAtomIdx()]
    start_pos = Point3D(start_pos[0], start_pos[1], start_pos[2])
    end_index = bond.GetEndAtomIdx()
    end_pos = pos[bond.GetEndAtomIdx()]
    end_pos = Point3D(end_pos[0], end_pos[1], end_pos[2])
    radius = 0.2
    color = Color(120, 120, 120)
    index = bond.GetIdx()
    edge = Edge(index, start_index, start_pos, end_index, end_pos, radius, color)
    scene.add_edge(edge)

svg_str = scene.to_svg()

with open("./out/test.svg", "w") as f:
    f.write(svg_str)