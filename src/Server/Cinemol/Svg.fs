namespace Cinemol

module Svg =

    open Cinemol.Types

    let header ((xMin, yMin, width, height): ViewBox) =
        $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        \n<svg\
        \n\tid=\"Layer_1\"\
        \n\txmlns=\"http://www.w3.org/2000/svg\"\
        \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
        \n>"

    let writeAtomStyle ((idx, _, _): AtomInfo) =  $"\n.atom-{idx}{{fill:url(#radial-gradient-{idx});}}"

    let writeAtomDefs pov unitSize ((idx, atom, coords): AtomInfo) =
        let distZ = abs (pov.Z - coords.Z)
        let ratio = unitSize / distZ
        let radius = atom.Radius * ratio
        let d1, d2, d3, d4, d5 = atom.Color.Gradient
        $"\n<radialGradient\
        \n\tid=\"radial-gradient-{idx}\"\
        \n\tcx=\"{Operators.string coords.X}\"\
        \n\tcy=\"{Operators.string coords.Y}\"\
        \n\tfx=\"{Operators.string coords.X}\"\
        \n\tfy=\"{Operators.string coords.Y}\"\
        \n\tr=\"{Operators.string (radius + 0.4)}\"\
        \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{Operators.string dr1}\" stop-color=\"rgb({Operators.string d1.R},{Operators.string d1.G},{Operators.string d1.B})\"/>\
        \n<stop offset=\"{Operators.string dr2}\" stop-color=\"rgb({Operators.string d2.R},{Operators.string d2.G},{Operators.string d2.B})\"/>\
        \n<stop offset=\"{Operators.string dr3}\" stop-color=\"rgb({Operators.string d3.R},{Operators.string d3.G},{Operators.string d3.B})\"/>\
        \n<stop offset=\"{Operators.string dr4}\" stop-color=\"rgb({Operators.string d4.R},{Operators.string d4.G},{Operators.string d4.B})\"/>\
        \n<stop offset=\"{Operators.string dr5}\" stop-color=\"rgb({Operators.string d5.R},{Operators.string d5.G},{Operators.string d5.B})\"/>\
        \n</radialGradient>"

    let writeAtom pov unitSize ((idx, atom, coords): AtomInfo) =
        let distZ = abs (pov.Z - coords.Z)
        let ratio = unitSize / distZ
        let radius = atom.Radius * ratio
        $"\n<circle\
        \n\tclass=\"atom-{idx}\"\
        \n\tcx=\"{Operators.string coords.X}\"\
        \n\tcy=\"{Operators.string coords.Y}\"\
        \n\tr=\"{Operators.string radius}\"\
        \n/>"

    let writeAtomArc pov unitSize ((idx, atom, coords): AtomInfo) (start: Coords) (final: Coords) : string =
        let distZ = abs (pov.Z - coords.Z)
        let ratio = unitSize / distZ
        let radius = atom.Radius * ratio
        let arcSweep = "1"
        $"<path\
        \n\tclass=\"atom-{idx}\"
        \n\td=\"M {start.X} {start.Y} A {radius} {radius} 0 {arcSweep} {0} {final.X} {final.Y} L {start.X} {start.Y}\"
        \n/>"

    let writeAtomsStyle (atoms: AtomInfo[]) : string =  atoms |> Array.map writeAtomStyle |> String.concat ""
    let writeAtomsDefs pov unitSize (atoms: AtomInfo[]) : string = atoms |> Array.map (writeAtomDefs pov unitSize) |> String.concat ""
    let writeAtoms pov unitSize (atoms: AtomInfo[]) : string = atoms |> Array.map (writeAtom pov unitSize) |> String.concat ""

    let add (s: string) (sb: System.Text.StringBuilder) : System.Text.StringBuilder = sb.Append(s)
    let stringify (sb: System.Text.StringBuilder) : string = sb.ToString()

    let writeSVG viewBox depiction pov unitSize atoms : string =
        System.Text.StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyle atoms)
        |> add "\n</style>"
        |> add (writeAtomsDefs pov unitSize atoms)
        |> add "\n</defs>"
        |> add (writeAtoms pov unitSize atoms)
        |> add "\n</svg>"
        |> stringify