#!/usr/bin/env python3
import typing as ty 
import math 
import dataclasses

# on draggable SVG: 
# https://www.petercollingridge.co.uk/tutorials/svg/interactive/dragging/

# on clipping:
#   circle decomposition: 
#   http://complexdan.com/svg-circleellipse-to-path-converter/
#   https://stackoverflow.com/questions/29864022/drawing-parts-of-circles-circumference-in-svg
#   https://stackoverflow.com/questions/5736398/how-to-calculate-the-svg-path-for-an-arc-of-a-circle/24569190#24569190
#   https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/pointer-events
#   https://developer.mozilla.org/en-US/docs/Web/SVG/Element/clipPath

def svg_header():
    return '<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg viewBox="-400 -400 800 800" version="1.1"\nxmlns="http://www.w3.org/2000/svg" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n'
def svg_trailer():
    return '</svg>'
def svg_circle(cx: float, cy: float, r: float, color: str) -> str: 
    return f'<circle cx=\"{cx}\" cy=\"{cy}\" r=\"{r}\" style=\"fill: {color};\"/>\n'

def polar_to_cartesian(cx: float, cy: float, r: float, angle_in_degrees: float) -> ty.Tuple[float, float]:
    angle_in_radians = (angle_in_degrees - 90) * math.pi / 180
    return (cx+(r*math.cos(angle_in_radians)), cy+(r*math.sin(angle_in_radians)))

def describe_arc(x: float, y: float, r: float, start_angle: float, end_angle: float, color: str) -> float:
    start_x, start_y = polar_to_cartesian(x, y, r, end_angle)
    end_x, end_y = polar_to_cartesian(x, y, r, start_angle)
    if (end_angle - start_angle) <= 180: arc_sweep = "0"
    else: arc_sweep = "1"
    # return f'<path d=\"M {start_x} {start_y} A {r} {r} 0 {arc_sweep} 0 {end_x} {end_y} L {x} {y} L {start_x} {start_y}\" style=\"fill: {color}; stroke: #446688; stroke-width: 0;\"/>'
    return f'<path d=\"M {start_x} {start_y} A {r} {r} 0 {arc_sweep} 0 {end_x} {end_y} Q 95 10 {start_x} {start_y}\" style=\"fill: {color}; stroke: #446688; stroke-width: 0;\"/>'
    

Radius = float 

@dataclasses.dataclass
class Point:
    x: float 
    y: float 
    z: float 

class Sphere(ty.Protocol):
    def __call__(self, p: Point) -> float: ... 

class Plane(ty.Protocol):
    def __call__(self, p: Point) -> float: ...

def sphere(c: Point, r: Radius) -> Sphere:
    return lambda p: (p.x-c.x)**2 + (p.y-c.y)**2 + (p.z-c.z)**2 - r**2 

def intersect_spheres(a: Sphere, b: Sphere) -> Plane:
    return lambda p: a(p) - b(p)

def main() -> None:
    sphere_a = sphere(Point(0, 0, 0), 1)
    sphere_b = sphere(Point(1, 0, 0), 1)
    
    # Plane will return 0.0 if Point is on plane
    plane = intersect_spheres(sphere_a, sphere_b)
    point = Point(0.5, 0, 0)
    print(plane(point))

    with open("test.svg", "w") as fo:
        fo.write(svg_header())
        fo.write(svg_circle(0, 0, 90, "#ff0000"))
        # fo.write(svg_circle(-60, 35, 57, "#ccc"))
        # fo.write(svg_circle(60, 35, 57, "#ccc"))
        fo.write(describe_arc(-60.5, 34.5, 57, 150, 330, "#ccc"))
        fo.write(describe_arc(60.5, 34.5, 57, 30, 210, "#ccc"))
        fo.write(svg_trailer())

# NOTE: what do we do with multiple intersections?



if __name__ == "__main__":
    main()
