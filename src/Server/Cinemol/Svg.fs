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

    let writeAtomStyle atom =  $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}"

    let writeAtomDefs atom =
        let d1, d2, d3, d4, d5 = atom.Type.Color.Gradient
        $"\n<radialGradient\
        \n\tid=\"radial-gradient-{atom.Index}\"\
        \n\tcx=\"{Operators.string atom.ProjectedCenter.X}\"\
        \n\tcy=\"{Operators.string atom.ProjectedCenter.Y}\"\
        \n\tfx=\"{Operators.string atom.ProjectedCenter.X}\"\
        \n\tfy=\"{Operators.string atom.ProjectedCenter.Y}\"\
        \n\tr=\"{Operators.string (atom.ProjectedRadius + 0.4)}\"\
        \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{Operators.string diffusionRate1}\" \
        stop-color=\"rgb({Operators.string d1.R},{Operators.string d1.G},{Operators.string d1.B})\"/>\
        \n<stop offset=\"{Operators.string diffusionRate2}\" \
        stop-color=\"rgb({Operators.string d2.R},{Operators.string d2.G},{Operators.string d2.B})\"/>\
        \n<stop offset=\"{Operators.string diffusionRate3}\" \
        stop-color=\"rgb({Operators.string d3.R},{Operators.string d3.G},{Operators.string d3.B})\"/>\
        \n<stop offset=\"{Operators.string diffusionRate4}\" \
        stop-color=\"rgb({Operators.string d4.R},{Operators.string d4.G},{Operators.string d4.B})\"/>\
        \n<stop offset=\"{Operators.string diffusionRate5}\" \
        stop-color=\"rgb({Operators.string d5.R},{Operators.string d5.G},{Operators.string d5.B})\"/>\
        \n</radialGradient>"

    let writeAtom (atom: AtomInfo) (otherAtoms: AtomInfo[]) =
        let clipping =
            otherAtoms
            |> Array.map (fun a -> atom.Intersect a)
            |> Array.filter (fun c -> match c with | NoIntersection | IntersectionPoint _ -> false | _ -> true)

        match clipping with
        // No clipping
        | xs when Array.length xs = 0 ->
            $"\n<circle\
            \n\tclass=\"atom-{atom.Index}\"\
            \n\tcx=\"{Operators.string atom.ProjectedCenter.X}\"\
            \n\tcy=\"{Operators.string atom.ProjectedCenter.Y}\"\
            \n\tr=\"{Operators.string atom.ProjectedRadius}\"\
            \n/>"

        // Clipping
        | xs when Array.forall (fun x -> match x with | Eclipsed -> false | _ -> true) xs ->
            ""  // TODO

        // Atom is eclipsed, will not be drawn
        | _ -> ""

    let writeAtomArc atom (start: Point) (final: Point) : string =
        $"<path\
        \n\tclass=\"atom-{atom.Index}\"
        \n\td=\"M {start.X} {start.Y} A {atom.ProjectedRadius} {atom.ProjectedRadius} 0 1 {0} {final.X} {final.Y} L {start.X} {start.Y}\"
        \n/>"

    let writeAtomsStyle (atoms: AtomInfo[]) : string =  atoms |> Array.map writeAtomStyle |> String.concat ""
    let writeAtomsDefs (atoms: AtomInfo[]) : string = atoms |> Array.map writeAtomDefs |> String.concat ""
    let writeAtoms (atoms: AtomInfo[]) : string =
        atoms
        |> Array.map (fun atom ->
            let otherAtoms = atoms |> Array.filter (fun a -> a.Index <> atom.Index)
            writeAtom atom otherAtoms)
        |> String.concat ""

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