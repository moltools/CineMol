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

// =====================================================================================================================
// Style
// =====================================================================================================================
let writeAtomDefs (atom: ProjectedAtomInfo) : string =
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

let writeAtomsDefs (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map writeAtomDefs
    |> String.concat ""

let writeAtomsStyle (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}")
    |> String.concat ""

let writeBondDefs (atoms: ProjectedAtomInfo[]) (bond: BondInfo) : string =
    let atoms = Array.toList atoms

    let rec findAtom (l: ProjectedAtomInfo list) atomIdx : ProjectedAtomInfo option =
        match l with
        | [] -> None
        | x::xs -> if x.Index = atomIdx then Some x else findAtom xs atomIdx

    match findAtom atoms bond.Start, findAtom atoms bond.End with
    | Some s, Some e ->
        let r3s, g3s, b3s = (getAtomColor CPK s.AtomType).Diffuse atomColorGradient.[2]
        let r3e, g3e, b3e = (getAtomColor CPK e.AtomType).Diffuse atomColorGradient.[2]
        $"\n<linearGradient\
        \n\tid=\"linear-gradient-{bond.Index}-atom-{s.Index}\"\
        \n\tx1=\"{floatToStr s.Center.X}\"\
        \n\ty1=\"{floatToStr s.Center.Y}\"\
        \n\tx2=\"{floatToStr e.Center.X}\"\
        \n\ty2=\"{floatToStr e.Center.Y}\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[2])}\" \
        stop-color=\"rgb({floatToStr r3s},{floatToStr g3s},{floatToStr b3s})\"/>\
        \n</linearGradient>
        \n<linearGradient\
        \n\tid=\"linear-gradient-{bond.Index}-atom-{e.Index}\"\
        \n\tx1=\"{floatToStr s.Center.X}\"\
        \n\ty1=\"{floatToStr s.Center.Y}\"\
        \n\tx2=\"{floatToStr e.Center.X}\"\
        \n\ty2=\"{floatToStr e.Center.Y}\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{floatToStr (1.0 - atomColorGradient.[2])}\" \
        stop-color=\"rgb({floatToStr r3e},{floatToStr g3e},{floatToStr b3e})\"/>\
        \n</linearGradient>"
    | _ -> ""

let writeBondsDefs (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    bonds
    |> Array.map (writeBondDefs atoms)
    |> String.concat ""

let writeBondsStyle (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    bonds
    |> Array.map (fun bond ->
        $"\n.bond-{bond.Index}-atom-{bond.Start}{{fill:url(#linear-gradient-{bond.Index}-atom-{bond.Start});}}
        \n.bond-{bond.Index}-atom-{bond.End}{{fill:url(#linear-gradient-{bond.Index}-atom-{bond.End});}}"
    )
    |> String.concat ""

// =====================================================================================================================
// Shapes
// =====================================================================================================================
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
        \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
        \n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} {writeArc p2 atom.Radius} {arcs} L {floatToStr p1.X} {floatToStr p1.Y}\"
        \n/>"
    | [] -> ""

let drawAtomFilled (atom: ProjectedAtomInfo) : string =
    $"\n<circle\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
    \n\tcx=\"{floatToStr atom.Center.X}\"\
    \n\tcy=\"{floatToStr atom.Center.Y}\"\
    \n\tr=\"{floatToStr atom.Radius}\"\
    \n/>"

let writeAtomFilled (atom: ProjectedAtomInfo) : string =
    match atom.Clipping with
    | [] -> drawAtomFilled atom
    | _ -> reconstructShape atom

let writeAtomsFilled (atoms: ProjectedAtomInfo[]) : string =
    atoms
    |> Array.map (fun atom -> writeAtomFilled atom)
    |> String.concat ""

let writeAtomsWire (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let bonds = Array.toList bonds

    let rec findAtom (l: ProjectedAtomInfo list) atomIdx : ProjectedAtomInfo option =
        match l with
        | [] -> None
        | x::xs -> if x.Index = atomIdx then Some x else findAtom xs atomIdx

    let findBonds (l: BondInfo list) atomIdx : BondInfo list =
        l |> List.filter (fun b -> b.Start = atomIdx)

    let mutable drawnAtoms = []

    let drawAtom (atom: ProjectedAtomInfo) : string =
        let r, g, b = (getAtomColor CPK atom.AtomType).RGB
        $"<circle\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\tstyle=\"fill:rgb({intToStr r},{intToStr g},{intToStr b})\"
        \n\tcx=\"{floatToStr atom.Center.X}\"\
        \n\tcy=\"{floatToStr atom.Center.Y}\"\
        \n\tr=\"{floatToStr 0.05}\"\
        \n/>"

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
                        | Some endAtom ->
                            let s, e = startAtom, endAtom
                            let m = s.Center.Midpoint e.Center
                            let sR, sG, sB = (getAtomColor CPK s.AtomType).RGB
                            let eR, eG, eB = (getAtomColor CPK e.AtomType).RGB
                            yield
                                $"<line\
                                \n\tx1=\"{s.Center.X}\"\
                                \n\tx2=\"{m.X}\"\
                                \n\ty1=\"{s.Center.Y}\"\
                                \n\ty2=\"{m.Y}\"\
                                \n\tstyle=\"stroke:rgb({intToStr sR},{intToStr sG},{intToStr sB});stroke-width:0.1\"/>"
                            yield
                                $"<line\
                                \n\tx1=\"{m.X}\"\
                                \n\tx2=\"{e.Center.X}\"\
                                \n\ty1=\"{m.Y}\"\
                                \n\ty2=\"{e.Center.Y}\"\
                                \n\tstyle=\"stroke:rgb({intToStr eR},{intToStr eG},{intToStr eB});stroke-width:0.1\"/>"
                        | None -> ()
    ]
    |> String.concat ""

type BondEnd = | Start | End

let writeBallAndStick (atoms: ProjectedAtomInfo[]) (bonds: BondInfo[]) : string =
    let atoms = Array.toList atoms
    let bonds = Array.toList bonds

    let drawAtom (atom: ProjectedAtomInfo) : string =
        $"<circle\
        \n\tclass=\"atom-{atom.Index}\"\
        \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
        \n\tcx=\"{floatToStr atom.Center.X}\"\
        \n\tcy=\"{floatToStr atom.Center.Y}\"\
        \n\tr=\"{floatToStr (atom.Radius / 2.0)}\"\
        \n/>"

    let drawSingleBond (bond: BondInfo) (s: ProjectedAtomInfo) (e: ProjectedAtomInfo) : string =
        let round = round 3
        let width = 0.15 * bond.Scaling
        let sProj: Point2D = { X = s.Center.X; Y = s.Center.Y }
        let eProj: Point2D = { X = e.Center.X; Y = e.Center.Y }

        let factor = 0.20
        let sProj: Point2D = {
            X = sProj.X + (factor * s.Radius) * (sProj.FindVector eProj).X
            Y = sProj.Y + (factor * s.Radius) * (sProj.FindVector eProj).Y
        }
        let eProj: Point2D = {
            X = eProj.X + (factor * e.Radius) * (eProj.FindVector sProj).X
            Y = eProj.Y + (factor * e.Radius) * (eProj.FindVector sProj).Y
        }

        [|
            Start, sProj, sProj.Midpoint eProj
            End, sProj.Midpoint eProj, eProj
        |]
        |> Array.map (fun (bondEnd, sProj, eProj) ->
            let slope = calcSlope sProj eProj

            // Construct cylinder for bond line
            let perpSlope = -1.0 * (1.0 / slope)  // Opposite reciprocal to get perpendicular slope
            let t =  width / Math.Sqrt (1.0 + Math.Pow(perpSlope, 2.0))
            let sTop: Point2D = { X = round (sProj.X + t); Y = round (sProj.Y + (perpSlope * t)) }
            let sBot: Point2D = { X = round (sProj.X - t); Y = round (sProj.Y - (perpSlope * t)) }
            let eTop: Point2D = { X = round (eProj.X + t); Y = round (eProj.Y + (perpSlope * t)) }
            let eBot: Point2D = { X = round (eProj.X - t); Y = round (eProj.Y - (perpSlope * t)) }
            let sSweep, eSweep = if sProj.Y > eProj.Y then 1, 0 else 0, 1

            let bondEnd =
                match bondEnd with
                | Start -> s.Index
                | End -> e.Index

            // Draw bond
            $"<path \
            \n\tclass=\"bond-{bond.Index}-atom-{bondEnd}\"
            \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
            \n\td= \
            \n\t\t\"M {floatToStr sTop.X} {floatToStr sTop.Y}\
            \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr sSweep} {floatToStr sBot.X} {floatToStr sBot.Y}\
            \n\t\tL {floatToStr eBot.X} {floatToStr eBot.Y}\
            \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr eSweep} {floatToStr eTop.X} {floatToStr eTop.Y}\
            \n\t\tL {floatToStr sTop.X} {floatToStr sTop.Y}\"\
            \n/>"
        )
        |> String.concat ""

    let drawDoubleBond (bond: BondInfo) (s: ProjectedAtomInfo) (e: ProjectedAtomInfo) : string =
        let round = round 3
        let sProj: Point2D = { X = s.Center.X; Y = s.Center.Y }
        let eProj: Point2D = { X = e.Center.X; Y = e.Center.Y }

        let factor = 0.20
        let sProj: Point2D = {
            X = sProj.X + (factor * s.Radius) * (sProj.FindVector eProj).X
            Y = sProj.Y + (factor * s.Radius) * (sProj.FindVector eProj).Y
        }
        let eProj: Point2D = {
            X = eProj.X + (factor * e.Radius) * (eProj.FindVector sProj).X
            Y = eProj.Y + (factor * e.Radius) * (eProj.FindVector sProj).Y
        }

        [|
            Start, sProj, sProj.Midpoint eProj
            End, sProj.Midpoint eProj, eProj
        |]
        |> Array.map (fun (bondEnd, sProj, eProj) ->
            let width = 0.10 * bond.Scaling
            let slope = calcSlope sProj eProj
            let perpSlope = -1.0 * (1.0 / slope)  // Opposite reciprocal to get perpendicular slope
            let t =  width / Math.Sqrt (1.0 + Math.Pow(perpSlope, 2.0))

            let aSideTop, aSideBot: Point2D * Point2D = (
                { X = round (sProj.X + t); Y = round (sProj.Y + (perpSlope * t)) },
                { X = round (eProj.X + t); Y = round (eProj.Y + (perpSlope * t)) }
            )
            let bSideTop, bSideBot: Point2D * Point2D = (
                { X = round (sProj.X - t); Y = round (sProj.Y - (perpSlope * t)) },
                { X = round (eProj.X - t); Y = round (eProj.Y - (perpSlope * t)) }
            )

            let constructCylinder p1 p2 =
                let width = 0.05 * bond.Scaling
                // Construct cylinder for bond line
                let slope = calcSlope p1 p2
                let perpSlope = -1.0 * (1.0 / slope)  // Opposite reciprocal to get perpendicular slope
                let t =  width / Math.Sqrt (1.0 + Math.Pow(perpSlope, 2.0))
                let sTop: Point2D = { X = round (p1.X + t); Y = round (p1.Y + (perpSlope * t)) }
                let sBot: Point2D = { X = round (p1.X - t); Y = round (p1.Y - (perpSlope * t)) }
                let eTop: Point2D = { X = round (p2.X + t); Y = round (p2.Y + (perpSlope * t)) }
                let eBot: Point2D = { X = round (p2.X - t); Y = round (p2.Y - (perpSlope * t)) }
                let sSweep, eSweep = if p1.Y > p2.Y then 1, 0 else 0, 1
                sTop, sBot, eTop, eBot, sSweep, eSweep

            let sTop1, sBot1, eTop1, eBot1, sSweep1, eSweep1 = constructCylinder aSideTop aSideBot
            let sTop2, sBot2, eTop2, eBot2, sSweep2, eSweep2 = constructCylinder bSideTop bSideBot

            let bondEnd =
                match bondEnd with
                | Start -> s.Index
                | End -> e.Index

            // Draw bond
            [
                $"<path \
                \n\tclass=\"bond-{bond.Index}-atom-{bondEnd}\"
                \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
                \n\td= \
                \n\t\t\"M {floatToStr sTop1.X} {floatToStr sTop1.Y}\
                \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr sSweep1} {floatToStr sBot1.X} {floatToStr sBot1.Y}\
                \n\t\tL {floatToStr eBot1.X} {floatToStr eBot1.Y}\
                \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr eSweep1} {floatToStr eTop1.X} {floatToStr eTop1.Y}\
                \n\t\tL {floatToStr sTop1.X} {floatToStr sTop1.Y}\"\
                \n/>"

                $"<path \
                \n\tclass=\"bond-{bond.Index}-atom-{bondEnd}\"
                \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"
                \n\td= \
                \n\t\t\"M {floatToStr sTop2.X} {floatToStr sTop2.Y}\
                \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr sSweep2} {floatToStr sBot2.X} {floatToStr sBot2.Y}\
                \n\t\tL {floatToStr eBot2.X} {floatToStr eBot2.Y}\
                \n\t\tA {floatToStr width} {floatToStr width} 0 0 {intToStr eSweep2} {floatToStr eTop2.X} {floatToStr eTop2.Y}\
                \n\t\tL {floatToStr sTop2.X} {floatToStr sTop2.Y}\"\
                \n/>"
            ]
            |> String.concat ""
        )
        |> String.concat ""

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
                        | Some endAtom ->
                            match atomBond.BondType with
                            | Single -> yield drawSingleBond atomBond startAtom endAtom
                            | Double -> yield drawDoubleBond atomBond startAtom endAtom
                            | _ -> ()
                        | None -> ()
    ]
    |> String.concat ""

// =====================================================================================================================
// Compile SVG
// =====================================================================================================================
let add (s: string) (sb: StringBuilder) : StringBuilder = sb.Append(s)

let stringify (sb: StringBuilder) : string = sb.ToString()

let writeSVG viewBox depiction (mol: ProjectedMolecule) : string =
    match depiction with
    | Filled ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyle mol.Atoms)
        |> add "\n</style>"
        |> add (writeAtomsDefs mol.Atoms)
        |> add "\n</defs>"
        |> add (writeAtomsFilled mol.Atoms)
        |> add "\n</svg>"
        |> stringify

    | BallAndStick ->
        StringBuilder()
        |> add (header viewBox)
        |> add "\n<defs>\n<style>"
        |> add (writeAtomsStyle mol.Atoms)
        |> add (writeBondsStyle mol.Atoms mol.Bonds)
        |> add "\n</style>"
        |> add (writeAtomsDefs mol.Atoms)
        |> add (writeBondsDefs mol.Atoms mol.Bonds)
        |> add "\n</defs>"
        |> add (writeBallAndStick mol.Atoms mol.Bonds)
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
