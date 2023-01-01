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

let writeAtomsDefs (atoms: ProjectedAtomInfo list) : string =
    atoms
    |> List.map writeAtomDefs
    |> String.concat ""

let writeAtomsStyle (atoms: ProjectedAtomInfo list) : string =
    atoms
    |> List.map (fun atom -> $"\n.atom-{atom.Index}{{fill:url(#radial-gradient-{atom.Index});}}")
    |> String.concat ""

let writeBondDefs (atoms: ProjectedAtomInfo list) (bond: BondInfo) : string =
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

let writeBondsDefs (atoms: ProjectedAtomInfo list) (bonds: BondInfo list) : string =
    bonds
    |> List.map (writeBondDefs atoms)
    |> String.concat ""

let writeBondsStyle (bonds: BondInfo list) : string =
    bonds
    |> List.map (fun bond ->
        $"\n.bond-{bond.Index}-atom-{bond.Start}{{fill:url(#linear-gradient-{bond.Index}-atom-{bond.Start});}}
        \n.bond-{bond.Index}-atom-{bond.End}{{fill:url(#linear-gradient-{bond.Index}-atom-{bond.End});}}"
    )
    |> String.concat ""

// =====================================================================================================================
// Shapes
// =====================================================================================================================
let reconstructShape (atom: ProjectedAtomInfo) : string =
    /// Modulo 2 pi
    let rec norm (rad: float) : float =
        let pi2 = 2.0 * Math.PI
        match rad with
        | x when x > pi2 -> norm (x - pi2)
        | x when x < 0.0 -> norm (x + pi2)
        | x -> x

    let dist ref p = p - ref |> norm

    /// Direction center to point on circle in radians
    let dir (p: Point2D) : float = Math.Atan2(p.Y - atom.Center.Y, p.X - atom.Center.X)

    /// Order clippings
    let mutable start: float option = None
    let orderedClipping =
        [ for { Line = (p1, p2) } in atom.Clipping do
            match start with
            | None ->
                let a1, a2 = dir p1, dir p2
                start <- Some a1
                (dist a1 a1, p1), (dist a1 a2, p2)
            | Some ref ->
                let a1, a2 = dir p1, dir p2
                /// Point on circle closest to ref comes first; flip around clipping if not the case
                if dist ref a1 < dist ref a2 then (dist ref a1, p1), (dist ref a2, p2)
                else (dist ref a2, p2), (dist ref a1, p1) ]
        |> List.sortBy (fun ((x, _), (_, _)) -> x)

    /// Merge clippings if they intersect
    let rec mergeClipping (merged: ((float * Point2D) * (float * Point2D)) list)
                          (cs: ((float * Point2D) * (float * Point2D)) list) :
                          ((float * Point2D) * (float * Point2D)) list =
        match cs with
        | c1::c2::rest ->
            let (a1,p1),(a2,p2) = c1
            let (a3,p3),(a4,p4) = c2
            if a2 > a4 then mergeClipping (merged @ [c1]) rest
            elif a2 < a3 then mergeClipping (merged @ [c1]) ([c2] @ rest)
            else
//                let i = intersectionBetweenLines (p1, p2) (p3, p4)
//                let c1New = (a1,p1),(a2,i)
//                let c2New = (a3,i),(a4,p4)
                let c1New = (a1,p1),(a2,p2)
                let c2New = (a2,p2),(a4,p4)
                mergeClipping (merged @ [c1New]) ([c2New] @ rest)
        | [c] -> merged @ [c]
        | [] -> merged

    let preppedClipping: Line list =
        mergeClipping [] orderedClipping
        |> List.map (fun ((_, p1),(_, p2)) -> p1, p2)

    let writeArc (final: Point2D) : string =
        let radius = atom.Radius
        let sweep = "1"
        $"A {floatToStr radius} {floatToStr radius} 0 {sweep} 0 {floatToStr final.X} {floatToStr final.Y}"

    /// Draw ordered and merged clipping paths
//    let mutable start: Point2D option = None
//    let rec writeSvg (s: string list) (cs: Line list) : string list =
//        match cs with
//        | [] -> s
//        | [(p1, p2)] when s.Length = 0 && cs.Length = 1 ->
//            let arc = writeArc p1 atom.Radius
//            s @ [$"\n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} L {floatToStr p2.X} {floatToStr p2.Y} {arc}\""]
//        | (p1, p2)::rest when s.Length = 0 ->
//            start <- Some p1
//            writeSvg [ $"\n\td=\"M {floatToStr p1.X} {floatToStr p1.Y} L {floatToStr p2.X} {floatToStr p2.Y} " ] rest
//        | (p1, p2)::rest ->
//            let connectToStart =
//                match start with
//                | Some pStart -> writeArc pStart atom.Radius
//                | None -> "" // will connect with straight line in this case ...
//            let sNew =
//                if rest.Length = 0 then
//                    s @ [
//                        writeArc p1 atom.Radius
//                        $" L {floatToStr p2.X} {floatToStr p2.Y} "
//                        connectToStart
//                        "\""
//                    ]
//                else s @ [ writeArc p1 atom.Radius; $"L {floatToStr p2.X} {floatToStr p2.Y} " ]
//            writeSvg sNew rest


    let cs = preppedClipping |> List.toArray
    let clipping =
        if cs.Length = 1 then
            let p1, p2 = cs.[0]
            [ $"M {floatToStr p1.X} {floatToStr p1.Y} L {floatToStr p2.X} {floatToStr p2.Y} {writeArc p1}" ]
        elif cs.Length = 2 then
            let p1, p2 = cs.[0]
            let p3, p4 = cs.[1]
            [ $"M {floatToStr p1.X} {floatToStr p1.Y} L {floatToStr p2.X} {floatToStr p2.Y} L {floatToStr p3.X} {floatToStr p3.Y} L {floatToStr p4.X} {floatToStr p4.Y} L {floatToStr p1.X} {floatToStr p1.Y} " ]
        else
            let p1, p2 = cs.[0]
            let p3, p4 = cs.[cs.Length - 1]
            [ $"M {floatToStr p1.X} {floatToStr p1.Y} L {floatToStr p2.X} {floatToStr p2.Y} " ]
            @ [ for a, b in cs.[1..cs.Length - 1] do $"L {floatToStr a.X} {floatToStr a.Y} L {floatToStr b.X} {floatToStr b.Y}" ]
            @ [ $"L {floatToStr p3.X} {floatToStr p3.Y} L {floatToStr p4.X} {floatToStr p4.Y} L {floatToStr p1.X} {floatToStr p1.Y} " ]

    [ $"\n<path\
    \n\tclass=\"atom-{atom.Index}\"\
    \n\tstyle=\"stroke:rgb(0,0,0);stroke-width:0.025\"" ]
    @ ["\n\td=\""]
    @ clipping
    @ ["\"/>"]
    |> String.concat ""


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

let writeAtomsFilled (atoms: ProjectedAtomInfo list) : string =
    atoms
    |> List.map (fun atom -> writeAtomFilled atom)
    |> String.concat ""

let writeAtomsWire (atoms: ProjectedAtomInfo list) (bonds: BondInfo list) : string =
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

let writeBallAndStick (atoms: ProjectedAtomInfo list) (bonds: BondInfo list) : string =
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
        |> add (writeBondsStyle mol.Bonds)
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
