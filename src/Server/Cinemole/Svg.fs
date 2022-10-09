namespace Cinemole

module Svg =

    open System

    open Cinemole.Helpers
    open Cinemole.Types
    open Cinemole.Geometry

    let header ((xMin, yMin, width, height): ViewBox) =
        $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        \n<svg\
        \n\tid=\"Layer_1\"\
        \n\txmlns=\"http://www.w3.org/2000/svg\"\
        \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
        \n>"

    let writeAtomStyle atom =  $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}"

    let writeAtomDefs atom =
        let d1, d2, d3, d4, d5 = atom.Type.Color.Gradient
        $"\n<radialGradient\
        \n\tid=\"radial-gradient-{atom.Index}\"\
        \n\tcx=\"{floatToStr atom.ProjectedCenter.X}\"\
        \n\tcy=\"{floatToStr atom.ProjectedCenter.Y}\"\
        \n\tfx=\"{floatToStr atom.ProjectedCenter.X}\"\
        \n\tfy=\"{floatToStr atom.ProjectedCenter.Y}\"\
        \n\tr=\"{floatToStr (atom.ProjectedRadius + 0.4)}\"\
        \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{floatToStr diffusionRate1}\" \
        stop-color=\"rgb({intToStr d1.R},{intToStr d1.G},{intToStr d1.B})\"/>\
        \n<stop offset=\"{floatToStr diffusionRate2}\" \
        stop-color=\"rgb({intToStr d2.R},{intToStr d2.G},{intToStr d2.B})\"/>\
        \n<stop offset=\"{floatToStr diffusionRate3}\" \
        stop-color=\"rgb({intToStr d3.R},{intToStr d3.G},{intToStr d3.B})\"/>\
        \n<stop offset=\"{floatToStr diffusionRate4}\" \
        stop-color=\"rgb({intToStr d4.R},{intToStr d4.G},{intToStr d4.B})\"/>\
        \n<stop offset=\"{floatToStr diffusionRate5}\" \
        stop-color=\"rgb({intToStr d5.R},{intToStr d5.G},{intToStr d5.B})\"/>\
        \n</radialGradient>"

    let alphaIsoscelesTriangle (p1: Point) (p2: Point) (r: Radius) : float =
        Math.Acos((2.0 * r ** 2.0 - (p1.Distance p2)) / (2.0 * r ** 2.0))

    // Sort edges based on how close they are together end-to-end based on respective alpha isoscele angle
    let sortEdges (r: Radius) (edges: (Point * Point)[]) : (Point * Point) list =
        // Orient edges so that first point in every edge is closer to the previous edge
        let rec orient (prev: (Point * Point) option) acc es : (Point * Point) list =
            match es with
            | (p1, p2)::tail ->
                match prev with
                | None -> orient (Some (p1, p2)) [(p1, p2)] tail
                | Some (prev1, prev2) ->
                    let p2prev1, p1prev2, p2prev2 = p2.Distance prev1, p1.Distance prev2, p2.Distance prev2
                    if p2prev2 < p2prev1 && p2prev2 < p1prev2 then orient (Some (p2, p1)) (acc @ [(p2, p1)]) tail
                    else orient (Some (p1, p2)) (acc @ [(p1, p2)]) tail
            | [] -> acc

        // Sort order of edges based on vicinity first edge
        match Array.toList edges with
        | (p1, p2)::tail ->
            if List.length tail = 0 then [(p1, p2)]
            elif List.length tail = 1 then [(p1, p2)] @ tail
            else
                tail
                |> List.map (fun (p3, p4) ->
                    let dist = (p1.Centroid p2).Distance (p3.Centroid p4)
                    p3, p4, dist)
                |> List.sortBy (fun (_, _, d) -> d)
                |> List.map (fun (p3, p4, _) -> (p3, p4))
        | [] -> []
        |> orient None []

    let writeArc (final: Point) (radius: Radius) : string =
        $"A {floatToStr radius} {floatToStr radius} 0 1 0 {floatToStr final.X} {floatToStr final.Y}"

    let reconstructShape (index: Index) (radius: Radius) (edges: (Point * Point)[]) : string =
        match Array.toList edges with
        | (p1, p2)::tail ->
            let arcs = [ for p3, p4 in tail do $" L {floatToStr p3.X} {floatToStr p3.Y} {writeArc p4 radius} " ] |> String.concat ""
            $"\n<path\
            \n\tclass=\"atom-{index}\"\
            \n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} {writeArc p2 radius} {arcs} L {floatToStr p1.X} {floatToStr p1.Y}\"
            \n/>"
        | [] -> ""

    let writeAtom (atom: AtomInfo) (atoms: AtomInfo[]) =
        // We only clip atom representation if it clips with atoms that are drawn earlier
        let otherAtoms = [|
            let mutable included = true
            for a in atoms do
                if a.Index = atom.Index then included <- false
                if included then yield a |]

        let clipping =
            otherAtoms
            |> Array.map (fun a -> atom.Intersect a)
            |> Array.filter (fun c -> match c with | NoIntersection | IntersectionPoint _ -> false | _ -> true)

        match clipping with
        // No clipping
        | xs when Array.length xs = 0 ->
            $"\n<circle\
            \n\tclass=\"atom-{atom.Index}\"\
            \n\tcx=\"{floatToStr atom.ProjectedCenter.X}\"\
            \n\tcy=\"{floatToStr atom.ProjectedCenter.Y}\"\
            \n\tr=\"{floatToStr atom.ProjectedRadius}\"\
            \n/>"

        // Clipping
        | xs when Array.forall (fun x -> match x with | Eclipsed -> false | _ -> true) xs ->
            [|
                for x in xs do
                    match x with
                    | IntersectionCircle (c, r, _) ->
                        match intersectionCircles atom.ProjectedCenter atom.ProjectedRadius c r with
                        | None -> ()
                        | Some (p1, p2) -> yield p1, p2
                    | _ -> ()
            |]
            |> reconstructShape atom.Index atom.ProjectedRadius

        // Eclipsed atom will not be drawn
        | _ -> ""

    let writeAtomsStyle (atoms: AtomInfo[]) : string =  atoms |> Array.map writeAtomStyle |> String.concat ""
    let writeAtomsDefs (atoms: AtomInfo[]) : string = atoms |> Array.map writeAtomDefs |> String.concat ""
    let writeAtoms (atoms: AtomInfo[]) : string = atoms |> Array.map (fun atom -> writeAtom atom atoms) |> String.concat ""

    let add (s: string) (sb: System.Text.StringBuilder) : System.Text.StringBuilder = sb.Append(s)
    let stringify (sb: System.Text.StringBuilder) : string = sb.ToString()

    let writeSVG viewBox depiction atoms : string =
        System.Text.StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyle atoms)
        |> add "\n</style>"
        |> add (writeAtomsDefs atoms)
        |> add "\n</defs>"
        |> add (writeAtoms atoms)
        |> add "\n</svg>"
        |> stringify