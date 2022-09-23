namespace Cinemol

module Geometry =

    open System

    open Cinemol.Types

    let abs (v : float) : float =
        (v ** 2.0) ** 0.5

    let dist (a : Coords) (b : Coords) : float =
        (Coords.Sum(Coords.Pow (a - b) 2.0)) ** 0.5

    let rotateAxisX (c : Coords) (rad : float option) : Coords =
        match rad with
        | Some r ->
            { X = c.X
              Y = c.Y * Math.Cos(r) - c.Z * Math.Sin(r)
              Z = c.Y * Math.Sin(r) + c.Z * Math.Cos(r) }
        | None -> c

    let rotateAxisY (c : Coords) (rad : float option) : Coords =
        match rad with
        | Some r ->
            { X = c.X * Math.Cos(r) + c.Z * Math.Sin(r)
              Y = c.Y
              Z = c.Z * Math.Cos(r) - c.X * Math.Sin(r) }
        | None -> c

    let rotateAxisZ (c : Coords) (rad : float option) : Coords =
        match rad with
        | Some r ->
            { X = c.X * Math.Cos(r) - c.Y * Math.Sin(r)
              Y = c.X * Math.Sin(r) + c.Y * Math.Cos(r)
              Z = c.Z }
        | None -> c

    let intersection_circles (c_p1: Coords) (r_p1: float) (c_p2: Coords) (r_p2: float) : (Coords * Coords) option =
        let d = Math.Sqrt((c_p2.X - c_p1.X) ** 2.0 + (c_p2.Y - c_p1.Y) ** 2.0)
        if d > (r_p1 + r_p2) then None  // non-intersecting
        elif d < abs (r_p1 - r_p1) then None  // one circle within other circle
        elif d = 0.0 && r_p1 = r_p2 then None  // coincident circles
        else
            let a = (r_p1 ** 2.0 - r_p2 ** 2.0 + d ** 2.0) / (2.0 * d)
            let h = Math.Sqrt(r_p1 ** 2.0 - a ** 2.0)
            let x2 = c_p1.X + a * (c_p2.X - c_p1.X) / d
            let y2 = c_p1.Y + a * (c_p2.Y - c_p1.Y) / d
            let x3 = x2 + h * (c_p2.Y - c_p1.Y) / d
            let y3 = y2 - h * (c_p2.X - c_p1.X) / d
            let x4 = x2 - h * (c_p2.Y - c_p1.Y) / d
            let y4 = y2 + h * (c_p2.X - c_p1.X) / d
            Some ({ X = x3; Y = y3; Z = 0.0 }, { X = x4; Y = y4; Z = 0.0 })