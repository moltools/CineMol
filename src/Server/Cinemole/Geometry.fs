namespace Cinemole

module Geometry =

    open System

    open Cinemole.Types

    let intersectionCircles (c_p1: Point) (r_p1: float) (c_p2: Point) (r_p2: float) : (Point * Point) option =
        let d = Math.Sqrt((c_p2.X - c_p1.X) ** 2.0 + (c_p2.Y - c_p1.Y) ** 2.0)
        if d > (r_p1 + r_p2) then None  // Non-intersecting
        elif d < abs (r_p1 - r_p1) then None  // One circle within other circle
        elif d = 0.0 && r_p1 = r_p2 then None  // Coincident circles
        else
            let d = d + 1E-05  // Make sure d is non-zero
            let a = (r_p1 ** 2.0 - r_p2 ** 2.0 + d ** 2.0) / (2.0 * d)
            let h = Math.Sqrt(r_p1 ** 2.0 - a ** 2.0)
            let x2 = c_p1.X + a * (c_p2.X - c_p1.X) / d
            let y2 = c_p1.Y + a * (c_p2.Y - c_p1.Y) / d
            let x3 = x2 + h * (c_p2.Y - c_p1.Y) / d
            let y3 = y2 - h * (c_p2.X - c_p1.X) / d
            let x4 = x2 - h * (c_p2.Y - c_p1.Y) / d
            let y4 = y2 + h * (c_p2.X - c_p1.X) / d
            Some ({ X = x3; Y = y3; Z = 0.0 }, { X = x4; Y = y4; Z = 0.0 })