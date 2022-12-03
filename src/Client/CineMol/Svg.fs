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

let writeAtom (atom: AtomInfo) (atoms: AtomInfo[]) =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.ProjectedCenter.X}\"\
    \n\tcy=\"{floatToStr atom.ProjectedCenter.Y}\"\
    \n\tr=\"{floatToStr atom.ProjectedRadius}\"\
    \n/>"

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

let writeSVG viewBox depiction atoms : string =
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
