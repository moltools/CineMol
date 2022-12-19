module Client.CineMol.Svg

open System
open System.Text

open Helpers
open Styles
open Types
open Geometry

let header ((xMin, yMin, width, height): ViewBox) =
    $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
    \n<svg\
    \n\tid=\"Layer_1\"\
    \n\txmlns=\"http://www.w3.org/2000/svg\"\
    \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
    \n>"

let writeAtomStyle atom =
    $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}"

let writeAtomDefs atom =
    let r1, g1, b1 = (getAtomColor CPK atom.AtomType).Diffuse atomColorGradient.[0]
    let r2, g2, b2 = (getAtomColor CPK atom.AtomType).Diffuse atomColorGradient.[1]
    let r3, g3, b3 = (getAtomColor CPK atom.AtomType).Diffuse atomColorGradient.[2]
    let r4, g4, b4 = (getAtomColor CPK atom.AtomType).Diffuse atomColorGradient.[3]
    let r5, g5, b5 = (getAtomColor CPK atom.AtomType).Diffuse atomColorGradient.[4]
    $"\n<radialGradient\
    \n\tid=\"radial-gradient-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.Center.X}\"\
    \n\tcy=\"{floatToStr atom.Center.Y}\"\
    \n\tfx=\"{floatToStr atom.Center.X}\"\
    \n\tfy=\"{floatToStr atom.Center.Y}\"\
    \n\tr=\"{floatToStr (atom.Radius + 0.4)}\"\
    \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
    \n\tgradientUnits=\"userSpaceOnUse\"\
    \n>\
    \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[0])}\" \
    stop-color=\"rgb({floatToStr r1},{floatToStr g1},{floatToStr b1})\"/>\
    \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[1])}\" \
    stop-color=\"rgb({floatToStr r2},{floatToStr g2},{floatToStr b2})\"/>\
    \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[2])}\" \
    stop-color=\"rgb({floatToStr r3},{floatToStr g3},{floatToStr b3})\"/>\
    \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[3])}\" \
    stop-color=\"rgb({floatToStr r4},{floatToStr g4},{floatToStr b4})\"/>\
    \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[4])}\" \
    stop-color=\"rgb({floatToStr r5},{floatToStr g5},{floatToStr b5})\"/>\
    \n</radialGradient>"

let alphaIsoscelesTriangle (p1: Point2D) (p2: Point2D) (r: Radius) : float =
    Math.Acos((2.0 * r ** 2.0 - (p1.Distance p2)) / (2.0 * r ** 2.0))

let writeArc (final: Point2D) (radius: Radius) : string =
    $"A {floatToStr radius} {floatToStr radius} 0 1 0 {floatToStr final.X} {floatToStr final.Y}"

let reconstructShape (atom: AtomInfo) (edges: (Point2D * Point2D)[]) : string =
    match Array.toList edges with
    | (p1, p2)::tail ->
        let arcs =
            [ for p3, p4 in tail do $" L {floatToStr p3.X} {floatToStr p3.Y} {writeArc p4 atom.Radius} " ]
            |> String.concat ""

        $"\n<path\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} {writeArc p2 atom.Radius} {arcs} L {floatToStr p1.X} {floatToStr p1.Y}\"
        \n/>"

    | [] -> ""

let drawAtom (atom: AtomInfo) : string =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.Center.X}\"\
    \n\tcy=\"{floatToStr atom.Center.Y}\"\
    \n\tr=\"{floatToStr atom.Radius}\"\
    \n/>"

let writeAtom (atom: AtomInfo) (otherAtoms: AtomInfo[]) =
    let clippingAtoms = [|
        for otherAtom in otherAtoms do
            match atom.Intersect otherAtom with
            | IntersectionCircle (p, r, v) -> yield (p, r, v)
            | _ -> () |]

    match clippingAtoms with
    // No clipping
    | cs when cs.Length = 0 ->
        drawAtom atom

    // Clipping
    | cs ->
            let intersections =
                cs
                |> Array.map (fun (p, r, v) ->
                    let center2D = { X = atom.Center.X; Y = atom.Center.Y }
                    let clipWithCenter2D = { X = p.X; Y = p.Y }
                    intersectionBetweenCircles center2D atom.Radius clipWithCenter2D r)
            let intersectionPoints = [| for intersection in intersections do match intersection with | None -> () | Some (p1, p2) -> yield (p1, p2) |]
            if intersectionPoints.Length = 0 then drawAtom atom
            else reconstructShape atom intersectionPoints

let writeAtomsStyle (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map writeAtomStyle
    |> String.concat ""

let writeAtomsDefs (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map writeAtomDefs
    |> String.concat ""

let writeAtoms (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> writeAtom atom atoms)
    |> String.concat ""

let add (s: string) (sb: StringBuilder) : StringBuilder = sb.Append(s)

let stringify (sb: StringBuilder) : string = sb.ToString()

let writeSVG viewBox depiction atoms pov : string =
    StringBuilder()
    |> add (header viewBox)
    |> add "\n<defs>\n<style>"
    |> add (writeAtomsStyle atoms)
    |> add "\n</style>"
    |> add (writeAtomsDefs atoms)
    |> add "\n</defs>"
    |> add (writeAtoms atoms)
    |> add "\n</svg>"
    |> stringify
