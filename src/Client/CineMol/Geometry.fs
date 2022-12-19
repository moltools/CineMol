module Client.CineMol.Geometry

open System

open Types

let intersectionBetweenCircles
    (c_p1: Point2D)
    (r_p1: float)
    (c_p2: Point2D)
    (r_p2: float)
    : (Point2D * Point2D) option =
    // Calculate the distance between the center of the two circles
    let d = Math.Sqrt((c_p2.X - c_p1.X) ** 2.0 + (c_p2.Y - c_p1.Y) ** 2.0)
    // Non-intersecting circles
    if d > (r_p1 + r_p2) then None
    // Coincident circles
    elif d = 0.0 && r_p1 = r_p2 then None
    // One circle within other circle
    elif d < abs (r_p1 - r_p1) || d < 1E-5 then None
    // Two intersection points
    else
        let a = (r_p1 ** 2.0 - r_p2 ** 2.0 + d ** 2.0) / (2.0 * d)
        let h1 = r_p1 ** 2.0 - a ** 2.0
        if h1 <= 0.0 then None
        else
            let h2 = Math.Sqrt(h1)
            let x2 = c_p1.X + a * (c_p2.X - c_p1.X) / d
            let y2 = c_p1.Y + a * (c_p2.Y - c_p1.Y) / d
            let x3 = x2 + h2 * (c_p2.Y - c_p1.Y) / d
            let y3 = y2 - h2 * (c_p2.X - c_p1.X) / d
            let x4 = x2 - h2 * (c_p2.Y - c_p1.Y) / d
            let y4 = y2 + h2 * (c_p2.X - c_p1.X) / d
            Some ({ X = x3; Y = y3 }, { X = x4; Y = y4 })

let clip (projAtom: AtomInfo) (projMol: Molecule) (atom: AtomInfo) (mol: Molecule) : Clipping list =
    let clippingAtoms = [|
        for projOtherAtom, otherAtom in Array.zip projMol.Atoms mol.Atoms do
            match atom.Intersect otherAtom with
            | IntersectionCircle _ ->
                match projAtom.Intersect projOtherAtom with
                | IntersectionCircle (p, r, _) -> yield (p, r)
                | _ -> ()
            | _ -> () |]

    match clippingAtoms with
    // No clipping
    | cs when cs.Length = 0 -> []
    // Clipping
    | cs ->
            let intersections =
                cs
                |> Array.map (fun (p, r) ->
                    let center2D = { X = projAtom.Center.X; Y = projAtom.Center.Y }
                    let clipWithCenter2D = { X = p.X; Y = p.Y }
                    intersectionBetweenCircles center2D projAtom.Radius clipWithCenter2D r)
            [ for intersection in intersections do
                match intersection with
                | None -> ()
                | Some (p1, p2) -> yield { Line = p1, p2; Sweep = One } ]
