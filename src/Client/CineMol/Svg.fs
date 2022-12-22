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

let writeAtomDefsFilled (atom: ProjectedAtomInfo) : string =
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

let reconstructShape (atom: ProjectedAtomInfo) : string =
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
        \n/>\n"
    | [] -> ""

let drawAtomFilled (atom: ProjectedAtomInfo) : string =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tcx=\"{floatToStr atom.Center.X}\"\
    \n\tcy=\"{floatToStr atom.Center.Y}\"\
    \n\tr=\"{floatToStr atom.Radius}\"\
    \n/>\n"

let writeAtomFilled (atom: ProjectedAtomInfo) : string =
    match atom.Clipping with
    | [] -> drawAtomFilled atom
    | _ -> reconstructShape atom

let writeAtomsStyleFilled (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}")
    |> String.concat ""

let writeAtomsDefsFilled (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map writeAtomDefsFilled
    |> String.concat ""

let writeAtomsFilled (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> writeAtomFilled atom)
    |> String.concat ""

let writeAtomsWire (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let rec findAtom (l: ProjectedAtomInfo list) idx =
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

let writeBallAndStick (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let bonds = Array.toList bonds

    let drawAtom (atom: ProjectedAtomInfo) : string =
        $"\n<circle\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\tcx=\"{floatToStr atom.Center.X}\"\
        \n\tcy=\"{floatToStr atom.Center.Y}\"\
        \n\tr=\"{floatToStr (atom.Radius / 2.0)}\"\
        \n/>\n"

    let drawBond (bond: BondInfo) (s: ProjectedAtomInfo) (e: ProjectedAtomInfo) : string =
        let width = 0.25 * bond.Scaling
        let sProj = { X = s.Center.X; Y = s.Center.Y }
        let eProj = { X = e.Center.X; Y = e.Center.Y }
        let slope = -1.0 * calcSlope sProj eProj
        let t =  width / Math.Sqrt (1.0 + Math.Pow(slope, 2.0))
        let sTop = { X = sProj.X + t; Y = sProj.Y + (slope * t) }
        let sBot = { X = sProj.X - t; Y = sProj.Y - (slope * t) }
        let eTop = { X = eProj.X + t; Y = eProj.Y + (slope * t) }
        let eBot = { X = eProj.X - t; Y = eProj.Y - (slope * t) }

        $"<path \
        \n\tclas=\"atom-{s.Index}\"
        \n\td= \
        \n\t\t\"M {floatToStr sTop.X} {floatToStr sTop.Y}\
        \n\t\tA {floatToStr width} {floatToStr width} 0 0 0 {floatToStr sBot.X} {floatToStr sBot.Y}\
        \n\t\tL {floatToStr eBot.X} {floatToStr eBot.Y}\
        \n\t\tA {floatToStr width} {floatToStr width} 0 0 0 {floatToStr eTop.X} {floatToStr eTop.Y}\
        \n\t\tL {floatToStr sTop.X} {floatToStr sTop.Y}\"\
        />\n"

    let rec findAtom (l: ProjectedAtomInfo list) atomIdx : ProjectedAtomInfo option =
        match l with
        | [] -> None
        | x::xs -> if x.Index = atomIdx then Some x else findAtom xs atomIdx

    let findBonds (l: BondInfo list) atomIdx : BondInfo list =
        l |> List.filter (fun b -> b.Start = atomIdx)

    let mutable drawnAtoms = []

    [
        for startAtom in atoms do
            yield drawAtom startAtom
            drawnAtoms <- drawnAtoms @ [ startAtom.Index ]
            match findBonds bonds startAtom.Index with
            | [] -> ()
            | atomBonds ->
                for atomBond in atomBonds do
                    if not (List.contains atomBond.End drawnAtoms) then
                        match findAtom atoms atomBond.End with
                        | Some endAtom -> yield drawBond atomBond startAtom endAtom
                        | None -> ()
    ]
    |> String.concat ""

let writeTube (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let bonds = Array.toList bonds

    let drawBond (bond: BondInfo) (s: ProjectedAtomInfo) (e: ProjectedAtomInfo) : string =
        let width = 0.25 * bond.Scaling
        let sProj = { X = s.Center.X; Y = s.Center.Y }
        let eProj = { X = e.Center.X; Y = e.Center.Y }
        let slope = -1.0 * calcSlope sProj eProj
        let t =  width / Math.Sqrt (1.0 + Math.Pow(slope, 2.0))
        let sTop = { X = sProj.X + t; Y = sProj.Y + (slope * t) }
        let sBot = { X = sProj.X - t; Y = sProj.Y - (slope * t) }
        let eTop = { X = eProj.X + t; Y = eProj.Y + (slope * t) }
        let eBot = { X = eProj.X - t; Y = eProj.Y - (slope * t) }

        $"<path \
        \n\tclas=\"atom-{s.Index}\"
        \n\td= \
        \n\t\t\"M {floatToStr sTop.X} {floatToStr sTop.Y}\
        \n\t\tA {floatToStr width} {floatToStr width} 0 0 0 {floatToStr sBot.X} {floatToStr sBot.Y}\
        \n\t\tL {floatToStr eBot.X} {floatToStr eBot.Y}\
        \n\t\tA {floatToStr width} {floatToStr width} 0 0 0 {floatToStr eTop.X} {floatToStr eTop.Y}\
        \n\t\tL {floatToStr sTop.X} {floatToStr sTop.Y}\"\
        />\n"

    let rec findAtom (l: ProjectedAtomInfo list) atomIdx : ProjectedAtomInfo option =
        match l with
        | [] -> None
        | x::xs -> if x.Index = atomIdx then Some x else findAtom xs atomIdx

    let findBonds (l: BondInfo list) atomIdx : BondInfo list =
        l |> List.filter (fun b -> b.Start = atomIdx)

    let mutable drawnAtoms = []

    [
        for startAtom in atoms do
            drawnAtoms <- drawnAtoms @ [ startAtom.Index ]
            match findBonds bonds startAtom.Index with
            | [] -> ()
            | atomBonds ->
                for atomBond in atomBonds do
                    if not (List.contains atomBond.End drawnAtoms) then
                        match findAtom atoms atomBond.End with
                        | Some endAtom -> yield drawBond atomBond startAtom endAtom
                        | None -> ()
    ]
    |> String.concat ""

let add (s: string) (sb: StringBuilder) : StringBuilder = sb.Append(s)

let stringify (sb: StringBuilder) : string = sb.ToString()

let writeSVG viewBox depiction (mol: ProjectedMolecule) : string =
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

    | BallAndStick ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyleFilled mol.Atoms)
        |> add "\n</style>"
        |> add (writeAtomsDefsFilled mol.Atoms)
        |> add "\n</defs>"
        |> add (writeBallAndStick mol.Atoms mol.Bonds)
        |> add "\n</svg>"
        |> stringify

    | Tube ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add "\n</style>"
        |> add "\n</defs>"
        |> add (writeTube mol.Atoms mol.Bonds)
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
