#!/usr/bin/env python3
# svg3d :: https://prideout.net/blog/svg_wireframes/
# Single-file Python library for generating 3D wireframes in SVG format.
# Copyright (c) 2019 Philip Rideout
# Distributed under the MIT License, see bottom of file.

# https://prideout.net/blog/occlusion_sorting/
# https://en.wikipedia.org/wiki/Polygon_mesh


import math
import numpy as np
import pyrr
import svgwrite

create_ortho = pyrr.matrix44.create_orthogonal_projection
create_perspective = pyrr.matrix44.create_perspective_projection
create_lookat = pyrr.matrix44.create_look_at
quaternion = pyrr.quaternion

from typing import NamedTuple, Callable, Sequence


class Viewport(NamedTuple):
    minx: float = -0.5
    miny: float = -0.5
    width: float = 1.0
    height: float = 1.0

    @classmethod
    def from_aspect(cls, aspect_ratio: float):
        return cls(-aspect_ratio / 2.0, -0.5, aspect_ratio, 1.0)

    @classmethod
    def from_string(cls, string_to_parse):
        args = [float(f) for f in string_to_parse.split()]
        return cls(*args)


class Camera(NamedTuple):
    view: np.ndarray
    projection: np.ndarray


class Mesh(NamedTuple):
    faces: np.ndarray
    shader: Callable[[int, float], dict] = None
    style: dict = None
    circle_radius: float = 0


class Scene(NamedTuple):
    meshes: Sequence[Mesh]

    def add_mesh(self, mesh: Mesh):
        self.meshes.append(mesh)


class View(NamedTuple):
    camera: Camera
    scene: Scene
    viewport: Viewport = Viewport()


class Engine:
    def __init__(self, views, precision=5):
        self.views = views
        self.precision = precision

    def render(self, filename, size=(512, 512), viewBox="-0.5 -0.5 1.0 1.0", **extra):
        drawing = svgwrite.Drawing(filename, size, viewBox=viewBox, **extra)
        self.render_to_drawing(drawing)
        drawing.save()

    def render_to_drawing(self, drawing):
        for view in self.views:
            projection = np.dot(view.camera.view, view.camera.projection)

            clip_path = drawing.defs.add(drawing.clipPath())
            clip_min = view.viewport.minx, view.viewport.miny
            clip_size = view.viewport.width, view.viewport.height
            clip_path.add(drawing.rect(clip_min, clip_size))

            for mesh in view.scene.meshes:
                g = self._create_group(drawing, projection, view.viewport, mesh)
                g["clip-path"] = clip_path.get_funciri()
                drawing.add(g)

    def _create_group(self, drawing, projection, viewport, mesh):
        faces = mesh.faces
        shader = mesh.shader or (lambda face_index, winding: {})
        default_style = mesh.style or {}

        # Extend each point to a vec4, then transform to clip space.
        faces = np.dstack([faces, np.ones(faces.shape[:2])])
        faces = np.dot(faces, projection)

        # Reject trivially clipped polygons.
        xyz, w = faces[:, :, :3], faces[:, :, 3:]
        accepted = np.logical_and(np.greater(xyz, -w), np.less(xyz, +w))
        accepted = np.all(accepted, 2)  # vert is accepted if xyz are all inside
        accepted = np.any(accepted, 1)  # face is accepted if any vert is inside
        degenerate = np.less_equal(w, 0)[:, :, 0]  # vert is bad if its w <= 0
        degenerate = np.any(degenerate, 1)  # face is bad if any of its verts are bad
        accepted = np.logical_and(accepted, np.logical_not(degenerate))
        faces = np.compress(accepted, faces, axis=0)

        # Apply perspective transformation.
        xyz, w = faces[:, :, :3], faces[:, :, 3:]
        faces = xyz / w

        # Sort faces from back to front.
        face_indices = self._sort_back_to_front(faces)
        faces = faces[face_indices]

        # Apply viewport transform to X and Y.
        faces[:, :, 0:1] = (1.0 + faces[:, :, 0:1]) * viewport.width / 2
        faces[:, :, 1:2] = (1.0 - faces[:, :, 1:2]) * viewport.height / 2
        faces[:, :, 0:1] += viewport.minx
        faces[:, :, 1:2] += viewport.miny

        # Compute the winding direction of each polygon.
        windings = np.zeros(faces.shape[0])
        if faces.shape[1] >= 3:
            p0, p1, p2 = faces[:, 0, :], faces[:, 1, :], faces[:, 2, :]
            normals = np.cross(p2 - p0, p1 - p0)
            np.copyto(windings, normals[:, 2])

        group = drawing.g(**default_style)

        # Create circles.
        if mesh.circle_radius > 0:
            for face_index, face in enumerate(faces):
                style = shader(face_indices[face_index], 0)
                if style is None:
                    continue
                face = np.around(face[:, :2], self.precision)
                for pt in face:
                    group.add(drawing.circle(pt, mesh.circle_radius, **style))
            return group

        # Create polygons and lines.
        for face_index, face in enumerate(faces):
            style = shader(face_indices[face_index], windings[face_index])
            if style is None:
                continue
            face = np.around(face[:, :2], self.precision)
            if len(face) == 2:
                group.add(drawing.line(face[0], face[1], **style))
            else:
                group.add(drawing.polygon(face, **style))

        return group

    def _sort_back_to_front(self, faces):
        z_centroids = -np.sum(faces[:, :, 2], axis=1)
        for face_index in range(len(z_centroids)):
            z_centroids[face_index] /= len(faces[face_index])
        return np.argsort(z_centroids)


def icosahedron():
    """Construct a 20-sided polyhedron"""
    faces = [
        (0, 1, 2),
        (0, 2, 3),
        (0, 3, 4),
        (0, 4, 5),
        (0, 5, 1),
        (11, 7, 6),
        (11, 8, 7),
        (11, 9, 8),
        (11, 10, 9),
        (11, 6, 10),
        (1, 6, 2),
        (2, 7, 3),
        (3, 8, 4),
        (4, 9, 5),
        (5, 10, 1),
        (6, 7, 2),
        (7, 8, 3),
        (8, 9, 4),
        (9, 10, 5),
        (10, 6, 1),
    ]
    verts = [
        (0.000, 0.000, 1.000),
        (0.894, 0.000, 0.447),
        (0.276, 0.851, 0.447),
        (-0.724, 0.526, 0.447),
        (-0.724, -0.526, 0.447),
        (0.276, -0.851, 0.447),
        (0.724, 0.526, -0.447),
        (-0.276, 0.851, -0.447),
        (-0.894, 0.000, -0.447),
        (-0.276, -0.851, -0.447),
        (0.724, -0.526, -0.447),
        (0.000, 0.000, -1.000),
    ]
    return verts, faces

def subdivide(verts, faces):
    """
    Subdivide each triangle into four triangles, pushing verts to the unit sphere.
    """
    triangles = len(faces)
    for faceIndex in range(triangles):

        # Create three new verts at the midpoints of each edge:
        face = faces[faceIndex]
        a, b, c = np.float32([verts[vertIndex] for vertIndex in face])
        verts.append(pyrr.vector.normalize(a + b))
        verts.append(pyrr.vector.normalize(b + c))
        verts.append(pyrr.vector.normalize(a + c))

        # Split the current triangle into four smaller triangles:
        i = len(verts) - 3
        j, k = i + 1, i + 2
        faces.append((i, j, k))
        faces.append((face[0], i, k))
        faces.append((i, face[1], j))
        faces[faceIndex] = (k, j, face[2])

    return verts, faces


def rgb(r, g, b):
    r = max(0.0, min(r, 1.0))
    g = max(0.0, min(g, 1.0))
    b = max(0.0, min(b, 1.0))
    return svgwrite.utils.rgb(r * 255, g * 255, b * 255)


def sphere(view, input_color=[1, 1, 0], steps=1, position=[0, 0, 0], radius=1):
    # Sphere Shell
    verts, faces = icosahedron()
    for _ in range(steps): verts, faces = subdivide(verts, faces)
    verts, faces = np.float32(verts), np.int32(faces)
    faces = verts[faces] + position

    # # Sphere Shell
    # new_verts, new_faces = icosahedron()
    # for _ in range(steps): verts, newfaces = subdivide(new_verts, new_faces)
    # new_verts, new_faces = np.float32(new_verts), np.int32(new_faces)
    # new_faces = new_verts[new_faces] + [0, 1, 0]

    # faces = np.vstack([faces, new_faces])

    # Sphere Lighting
    ones = np.ones(faces.shape[:2] + (1,))
    eyespace_faces = np.dstack([faces, ones])
    eyespace_faces = np.dot(eyespace_faces, view)[:, :, :3]
    shininess = 100
    L = pyrr.vector.normalize(np.float32([20, 20, 50]))
    E = np.float32([0, 0, 1])
    H = pyrr.vector.normalize(L + E)

    def frontface_shader(face_index, winding):
        if winding < 0:
            return None
        face = eyespace_faces[face_index]
        p0, p1, p2 = face[0], face[1], face[2]
        N = pyrr.vector.normalize(pyrr.vector3.cross(p1 - p0, p2 - p0))
        df = max(0, np.dot(N, L))
        sf = pow(max(0, np.dot(N, H)), shininess)
        color = df * np.float32(input_color) + sf * np.float32([1, 1, 1])
        color = np.power(color, 1.0 / 2.2)
        return dict(
            fill=rgb(*color), 
            fill_opacity="1.0", 
            stroke=rgb(*color), 
            stroke_width="0.001"
        )

    return Mesh(faces, frontface_shader)


def tube(view, input_color, steps, start, end, radius):
    # Sphere Shell
    # verts, faces = icosahedron()
    # verts, faces = np.float32(verts), np.int32(faces)
    # faces = verts[faces]

    num_points = 6
    # points = []
    # for idx in range(num_points):
    #     point = [
    #         radius * math.cos((idx * 2 * math.pi) / num_points),
    #         radius * math.sin((idx * 2 * math.pi) / num_points),
    #         0
    #     ]
    #     points.append(point)

    p1 = np.array(start)
    p2 = np.array(end)
    v3 = p1 - p2
    v3 = v3 / np.linalg.norm(v3)
    e = np.array([0, 0, 0])
    e[np.argmin(np.abs(v3))] = 1
    v1 = np.cross(e, v3)
    v1 = v1 / np.linalg.norm(v3)
    v2 = np.cross(v3, v1)

    s = np.pi/3
    p = p1 + radius*(np.cos(s)*v1 + np.sin(s)*v2)

    points_a = []
    points_b = []
    for idx in range(num_points):
        phi = idx * 2 * np.pi / num_points
        point = p1 + (radius * (math.cos(phi)*v1 + math.sin(phi)*v2))
        points_a.append(point)
        points_b.append(point + p2)


    faces = np.array([np.array(points_b)])

    faces = np.vstack([faces, np.array([points_a])])

    print(faces.shape)

    # Sphere Lighting
    ones = np.ones(faces.shape[:2] + (1,))
    eyespace_faces = np.dstack([faces, ones])
    eyespace_faces = np.dot(eyespace_faces, view)[:, :, :3]
    shininess = 100
    L = pyrr.vector.normalize(np.float32([20, 20, 50]))
    E = np.float32([0, 0, 1])
    H = pyrr.vector.normalize(L + E)

    def frontface_shader(face_index, winding):
        if winding < 0:
            return None
        face = eyespace_faces[face_index]
        p0, p1, p2 = face[0], face[1], face[2]
        N = pyrr.vector.normalize(pyrr.vector3.cross(p1 - p0, p2 - p0))
        df = max(0, np.dot(N, L))
        sf = pow(max(0, np.dot(N, H)), shininess)
        color = df * np.float32([0, 1, 0]) + sf * np.float32([1, 1, 1])
        color = np.power(color, 1.0 / 2.2)
        return dict(
            fill=rgb(*color), 
            fill_opacity="1.0", 
            stroke=rgb(*color), 
            stroke_width="0.001"
        )

    return Mesh(faces, frontface_shader)

    
def main():
    fn = "./out/sphere_default.svg"

    projection = create_perspective(fovy=25, aspect=1, near=10, far=200)
    view = create_lookat(eye=[25, 20, 60], target=[0, 0, 0], up=[0, 1, 0])
    camera = Camera(view, projection)
    scene = Scene([])

    a = [0, -2, 0]
    b = [0, 2, 0]

    # mesh = sphere(view, input_color=[0, 1, 1], steps=3, position=a, radius=2)
    # scene.add_mesh(mesh)
    # mesh = sphere(view, input_color=[1, 0, 1], steps=3, position=b, radius=2)
    # scene.add_mesh(mesh)
    mesh = tube(view, input_color=[0, 0, 1], steps=3, start=a, end=b, radius=1)
    scene.add_mesh(mesh)

    Engine([View(camera, scene)]).render(fn)

    print(0)

    # create_octahedron_pair("./out/platonic_octahedron.svg") # hexahedron


if __name__ == "__main__":
    main()
