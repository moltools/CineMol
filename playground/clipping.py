#!/usr/bin/env python3
from __future__ import annotations
from ctypes import c_char_p
from dataclasses import dataclass
from typing import Tuple
import math 

# on draggable SVG: 
# https://www.petercollingridge.co.uk/tutorials/svg/interactive/dragging/

# on clipping:
#   circle decomposition: 
#   http://complexdan.com/svg-circleellipse-to-path-converter/
#   https://stackoverflow.com/questions/29864022/drawing-parts-of-circles-circumference-in-svg
#   https://stackoverflow.com/questions/5736398/how-to-calculate-the-svg-path-for-an-arc-of-a-circle/24569190#24569190
#   https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/pointer-events
#   https://developer.mozilla.org/en-US/docs/Web/SVG/Element/clipPath

Radius = float 
Color = str

@dataclass
class Point2D:
    x: float 
    y: float 

    def dist(s, o: Point2D) -> float:
        return math.sqrt((o.x - s.x) ** 2 + (o.y - s.y) ** 2)

@dataclass
class Point3D:
    x: float 
    y: float 
    z: float 

    def dist(s, o: Point3D) -> float:
        return math.sqrt((o.x - s.x) ** 2 + (o.y - s.y) ** 2 + (o.z - s.z) ** 2)

@dataclass
class Atom:
    type: str 
    center: Point3D

    def radius(self) -> Radius:
        if self.type == "H": return 0.8
        elif self.type == "O": return 1.0
        else: raise ValueError(f"unknown atom type: {self.type}")

    def color(self) -> Color:
        if self.type == "H": return "#ccc"
        elif self.type == "O": return "#ff0000"
        else: raise ValueError(f"unknown atom type: {self.type}")

def svg_header(w: float, h: float): return f'<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg viewBox="{-(w/2)} {-(h/2)} {w} {h}" version="1.1"\nxmlns="http://www.w3.org/2000/svg" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n'
def svg_trailer(): return '</svg>'
def svg_circle(c: Point2D, r: float, color: str) -> str:  return f'<circle cx=\"{c.x}\" cy=\"{c.y}\" r=\"{r}\" style=\"fill: {color};\"/>\n'

def circles_intersect(
    c_a: Point2D, 
    r_a: Radius, 
    c_b: Point2D, 
    r_b: Radius
) -> bool:
    d = math.sqrt((c_a.x - c_b.x) ** 2 + (c_a.y - c_b.y) ** 2)
    if d > r_a + r_b: return False # Non-intersecting 
    if d < abs(r_a - r_b): return False # One circle within other circle
    if d == 0 and r_a == r_b: return False # Coincident circles
    else: return True # Intersection

def intersection_circles(
    c_p1: Point2D, 
    r_p1: Radius, 
    c_p2: Point2D, 
    r_p2: Radius
) -> Tuple[Point2D, Point2D]:
    d = math.sqrt((c_p2.x - c_p1.x) ** 2 + (c_p2.y - c_p1.y) **2)
    a = (r_p1 ** 2 - r_p2 ** 2 + d ** 2) / (2 * d)
    h = math.sqrt(r_p1 ** 2 - a ** 2)
    x2 = c_p1.x + a * (c_p2.x - c_p1.x) / d   
    y2 = c_p1.y + a * (c_p2.y - c_p1.y) / d   
    x3 = x2 + h * (c_p2.y - c_p1.y) / d     
    y3 = y2 - h * (c_p2.x - c_p1.x) / d 
    x4 = x2 - h * (c_p2.y - c_p1.y) / d
    y4 = y2 + h * (c_p2.x - c_p1.x) / d
    return (Point2D(x3, y3), Point2D(x4, y4))

def polar_to_cartesian(c: Point2D, r: float, angle_in_degrees: float) -> Point2D:
    angle_in_radians = (angle_in_degrees - 90) * math.pi / 180
    return Point2D(c.x + (r*math.cos(angle_in_radians)), c.y + (r*math.sin(angle_in_radians)))

# def describe_arc(c: Point2D, r: Radius, start_angle: float, end_angle: float, color: str) -> str:
def describe_arc(c: Point2D, r: Radius, start: Point2D, end: Point2D, color: str) -> str:
    # start = polar_to_cartesian(c, r, end_angle)
    # end = polar_to_cartesian(c, r, start_angle)
    # if (end_angle - start_angle) <= 180: arc_sweep = "0"
    # else: arc_sweep = "1"
    arc_sweep = "1"
    return f'<path d=\"M {start.x} {start.y} A {r} {r} 0 {arc_sweep} 0 {end.x} {end.y} L {start.x} {start.y}\" style=\"fill: {color}; stroke: #446688; stroke-width: 0;\"/>'

def angle_triangle(c: Point2D, p1: Point2D, p2: Point2D) -> Tuple[float, float]:
    c_p1, c_p2, p1_p2 = c.dist(p1), c.dist(p2), p1.dist(p2)
    d = math.degrees(math.acos((c_p1 ** 2 + c_p2 ** 2 - p1_p2 ** 2) / (2 * c_p1 * c_p2)))
    return (d, 360 - d)

def main() -> None:
    with open("test.svg", "w") as fo:
        fo.write(svg_header(10, 10))

        # Draw oxygen
        O = Atom("O", Point3D(-0.0110, 0.9628, 0.0073))
        center2D_O = Point2D(O.center.x, O.center.y)
        fo.write(svg_circle(center2D_O, O.radius(), O.color()))

        # Drag hydrogens
        hydrogens = [Atom("H", Point3D(0.0021, -0.0041, 0.0020)),  Atom("H", Point3D(0.8669, 1.3681, 0.0011))]
        for H in hydrogens: 
            center2D_H = Point2D(H.center.x, H.center.y)
            if circles_intersect(center2D_O, O.radius(), center2D_H, H.radius()):
                p1, p2 = intersection_circles(center2D_O, O.radius(), center2D_H, H.radius())
                fo.write(describe_arc(center2D_H, H.radius(), p2, p1, H.color()))
            else:
                fo.write(svg_circle(center2D_H, H.radius(), H.color()))

        fo.write(svg_trailer())

if __name__ == "__main__":
    main()
