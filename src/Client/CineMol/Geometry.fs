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

let calcSlope (p1: Point2D) (p2: Point2D) : float =
    (p2.Y - p1.Y) / (p2.X - p1.X)

let sameSideOfLine (line: Point2D * Point2D) (p1: Point2D) (p2: Point2D) : bool =
    let l1, l2 = line
    let d (p: Point2D) = (p.X - l1.X) * (l2.Y - l1.Y) - (p.Y - l1.Y) * (l2.X - l1.X)
    match d p1 > 0.0, d p2 > 0.0 with
    | b1, b2 when b1 = b2 -> true
    | _ -> false

let clip (pov: Point3D) (persAtom: AtomInfo) (persMol: Molecule) (atom: AtomInfo) (mol: Molecule) : ClipPath list =
    // Only clip with atoms that are behind the atom of interest
    let distPovAtom = pov.Distance atom.Center
    let atomsForClipping =
        [ for a in mol.Atoms do
              if distPovAtom < pov.Distance a.Center then yield a ]
    let inds = atomsForClipping |> List.map (fun a -> a.Index)
    let persAtomsForClipping =
        [ for a in persMol.Atoms do
            if List.contains a.Index inds then yield a ]

    // Determine which atoms are actually clipping; store clipping line
    let clippingAtoms = [|
        for persOtherAtom, otherAtom in List.zip persAtomsForClipping atomsForClipping do
            match atom.Intersects otherAtom with
            | true ->
                match persAtom.Intersection persOtherAtom with
                | IntersectionCircle (p, r, _) ->
                    if pov.Distance p < distPovAtom then
                        // circle is more like ellips when not looking right
                        // towards the front, so we adjust the radius a bit
                        yield p, r * 0.8, persOtherAtom
                | _ -> ()
            | false -> ()
    |]

    match clippingAtoms with
    // No clipping
    | cs when cs.Length = 0 -> []

    // Clipping
    | cs ->
            let intersections =
                cs |> Array.map (fun (i_p, i_r, other) ->
                    let p1: Point2D = { X = persAtom.Center.X; Y = persAtom.Center.Y }
                    let p2: Point2D = { X = i_p.X; Y = i_p.Y }
                    intersectionBetweenCircles p1 persAtom.Radius p2 i_r,
                    p1, persAtom.Radius,
                    p2, i_r,
                    other, other.Radius)

            [ for intersection, p1, r1, p2, r2, other, rOther in intersections do
                match intersection with
                | None -> ()
                | Some (l1, l2) ->
                    let p =
                        match sameSideOfLine (l1, l2) p1 p2 with
                        | true ->
                            match p1.Distance p2 with
                            | x when x < r1 && r1 < rOther  -> IncludeBothSides
                            | x when x < r1 -> IncludeSide p1
                            | _ -> ExcludeSide p1
                        | false -> IncludeSide p1
                    yield { Line = l1, l2; SelectForSide = p } ]

let intersectionBetweenLines (l1: Point2D * Point2D) (l2: Point2D * Point2D) : Point2D option =
    let p1, p2 = l1
    let p3, p4 = l2
    let a1 = calcSlope p1 p2
    let a2 = calcSlope p3 p4
    if a1 = a2 then
        // lines are parallel
        None
    else
        let c1 = p1.Y - (a1 * p1.X)
        let c2 = p3.Y - (a2 * p3.X)
        if (c1 = infinity || c1 = -infinity) || (c2 = infinity || c2 = -infinity) then
            // Intersection happens somewhere off-screen
            None
        else
            let x = (c1 - c2) / (a2 - a1)
            let y = (a1 * x) + c1
            Some { X = x; Y = y }

