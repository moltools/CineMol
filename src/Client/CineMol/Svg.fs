module Client.CineMol.Svg

open System.Text

open Helpers
open Styles
open Types

let header ((xMin, yMin, width, height): ViewBox) =
    $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
    \n<svg\
    \n\tid=\"Layer_1\"\
    \n\txmlns=\"http://www.w3.org/2000/svg\"\
    \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
    \n>"

let writeAtomStyleFilled (atom: AtomInfo) : string =
    $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}"

let writeAtomDefsFilled (atom: AtomInfo) : string =
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

let reconstructShape (atom: AtomInfo) : string =
    let writeArc (final: Point2D) (radius: Radius) : string =
        let sweep = 1
        $"A {floatToStr radius} {floatToStr radius} 0 {sweep} 0 {floatToStr final.X} {floatToStr final.Y}"

    match atom.Clipping with
    | {Line = (p1, p2)}::tail ->
        let arcs =
            [ for {Line = (p3, p4)} in tail do $" L {floatToStr p3.X} {floatToStr p3.Y} {writeArc p4 atom.Radius} " ]
            |> String.concat ""
        $"\n<path\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} {writeArc p2 atom.Radius} {arcs} L {floatToStr p1.X} {floatToStr p1.Y}\"
        \n/>"
    | [] -> ""

let drawAtomFilled (atom: AtomInfo) : string =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.Center.X}\"\
    \n\tcy=\"{floatToStr atom.Center.Y}\"\
    \n\tr=\"{floatToStr atom.Radius}\"\
    \n/>"

let writeAtomFilled (atom: AtomInfo) : string =
    match atom.Clipping with
    | [] -> drawAtomFilled atom
    | _ -> reconstructShape atom

let writeAtomsStyleFilled (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map writeAtomStyleFilled
    |> String.concat ""

let writeAtomsDefsFilled (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map writeAtomDefsFilled
    |> String.concat ""

let writeAtomsFilled (atoms: AtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> writeAtomFilled atom)
    |> String.concat ""

let writeAtomsWire (atoms: AtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let rec findAtom (l: AtomInfo list) idx =
        match l with
        | [] -> None
        | x::xs -> if x.Index = idx then Some x else findAtom xs idx
    bonds
    |> Array.map (fun b ->
        match (findAtom atoms b.Start, findAtom atoms b.End) with
        | Some s, Some e ->
            $"<line\
            \n\tx1=\"{s.Center.X}\"\
            \n\tx2=\"{e.Center.X}\"\
            \n\ty1=\"{s.Center.Y}\"\
            \n\ty2=\"{e.Center.Y}\"\
            \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.1\"/>\n"
        | _ -> ""
    )
    |> String.concat ""

let add (s: string) (sb: StringBuilder) : StringBuilder = sb.Append(s)

let stringify (sb: StringBuilder) : string = sb.ToString()

let writeSVG viewBox depiction mol : string =
    match depiction with
    | Filled ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyleFilled mol.Atoms)
        |> add "\n</style>"
        |> add (writeAtomsDefsFilled mol.Atoms)
        |> add "\n</defs>"
        |> add (writeAtomsFilled mol.Atoms)
        |> add "\n</svg>"
        |> stringify

    | Wire ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add "\n</style>"
        |> add "\n</defs>"
        |> add (writeAtomsWire mol.Atoms mol.Bonds)
        |> add "\n</svg>"
        |> stringify

    | _ -> "fails_drawing"  // Other depictions not yet implemented
