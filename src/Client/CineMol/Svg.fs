module Client.CineMol.Svg

open System
open System.Text

open Helpers
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
    let d1, d2, d3, d4, d5 = atom.Type.Color.Gradient
    $"\n<radialGradient\
    \n\tid=\"radial-gradient-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.C.X}\"\
    \n\tcy=\"{floatToStr atom.C.Y}\"\
    \n\tfx=\"{floatToStr atom.C.X}\"\
    \n\tfy=\"{floatToStr atom.C.Y}\"\
    \n\tr=\"{floatToStr (atom.R + 0.4)}\"\
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

let writeArc (final: Point) (radius: Radius) : string =
    $"A {floatToStr radius} {floatToStr radius} 0 1 0 {floatToStr final.X} {floatToStr final.Y}"

let reconstructShape (atom: AtomInfo) (edges: (Point * Point)[]) : string =
    match Array.toList edges with
    | (p1, p2)::tail ->
        let arcs =
            [ for p3, p4 in tail do $" L {floatToStr p3.X} {floatToStr p3.Y} {writeArc p4 atom.R} " ]
            |> String.concat ""

        $"\n<path\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} {writeArc p2 atom.R} {arcs} L {floatToStr p1.X} {floatToStr p1.Y}\"
        \n/>"

    | [] -> ""

let drawAtom (atom: AtomInfo) : string =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.C.X}\"\
    \n\tcy=\"{floatToStr atom.C.Y}\"\
    \n\tr=\"{floatToStr atom.R}\"\
    \n/>"

let writeAtom (atom: AtomInfo) (otherAtoms: AtomInfo[]) =
    let clippingAtoms = [|
        for otherAtom in otherAtoms do
            match atom.Intersect otherAtom with
            | IntersectionCircle (p, r, v) -> yield (p, r, v)
            | _ -> () |]

    match clippingAtoms with

    /// No clipping. Draw atom as full circle.
    | cs when cs.Length = 0 -> drawAtom atom

    /// Clipping. Calculate clipped parts of atom before drawing.
    | cs ->
            let intersections = Array.map (fun (p, r, v) -> intersectionBetweenCircles atom.C atom.R p r) cs
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
