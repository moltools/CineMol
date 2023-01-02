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
let svgCircle (atom: ProjectedAtomInfo) : string =
    let str = floatToStr
    let x, y, r = str atom.Center.X, str atom.Center.Y, str atom.Radius
    $"<circle class=\"atom-{atom.Index}\" style=\"stroke:rgb(0,0,0);stroke-width:0.025\" cx=\"{x}\" cy=\"{y}\" r=\"{r}\"/>"

type Path = PathElement list

and Flag = | Zero | One
    with
    override x.ToString () : string =
        match x with | Zero -> "0" | One -> "1"

and PathElement =
    | Opening of className: string * style: string
    | Closing
    | Start of point: Point2D
    | Line of point: Point2D
    | Arc of radius: float * point: Point2D * largeArcFlag: Flag * sweepFlag: Flag
    with
    member x.ToSvg () : string =
        let str = floatToStr
        match x with
        | Opening (className, style) -> $"<path class=\"{className}\" style=\"{style}\" d=\""
        | Closing -> "\"/>"
        | Start p -> $"M {str p.X} {str p.Y}"
        | Line p -> $"L {str p.X} {str p.Y}"
        | Arc (r, p, largeArcFlag, sweepFlag) ->
            $"A {str r} {str r} 0 {largeArcFlag} {sweepFlag} {str p.X} {str p.Y}"

let pathToSvg (path: Path) : string =
    path |> List.map (fun x -> x.ToSvg()) |> String.concat " "

let drawAtom (atom: ProjectedAtomInfo) : string =
    let opening = Opening ($"atom-{atom.Index}", "stroke:rgb(0,0,0);stroke-width:0.025")
    let closing = Closing
    let start p = Start p
    let line p = Line p
    let arc p largeArcFlag sweepFlag = Arc (atom.Radius, p, largeArcFlag, sweepFlag)

    let determineFlags (p1: Point2D) (p2: Point2D) pos =
        match p1.Distance p2 < (2.0 * atom.Radius), pos with
        | true, SameSide -> Zero, One
        | true, OppositeSide -> One, One
        | false, SameSide -> One, Zero
        | false, OppositeSide -> Zero, Zero

    match atom.Clippings with
    /// No clipping; draw full circle for atom
    | [] -> svgCircle atom

    /// A single clipping; draw part atom circle with a chipped of clipping
    | [{Line = (p1, p2); AtomCentersPosition = pos}] ->
        let largeArcFlag, sweepFlag = determineFlags p1 p2 pos
        [ opening; start p1; line p2; arc p1 largeArcFlag sweepFlag; closing ]
        |> pathToSvg

    /// Multiple clippings; reconstruct shape by consecutively chipping of clippings
    | head::tail ->
        let { Line = (p1, p2); AtomCentersPosition = pos } = head
        let pi2 = 2.0 * Math.PI
        let rec normalize (rad: float) : float =
            match rad with
            | x when x > pi2 -> normalize (x - pi2)
            | x when x < 0.0 -> normalize (x + pi2)
            | x -> x

        let dist refTheta theta = theta - refTheta |> normalize

        /// Calculates direction center to point on circle (theta)
        let direction (p: Point2D) = Math.Atan2(p.Y - atom.Center.Y, p.X - atom.Center.X)

        /// Pick reference point
        let refTheta =
            let refTheta = direction p1
            if dist refTheta (direction p2) > pi2 then direction p2
            else refTheta

        /// Check for eclipsed clippings
        let sections =
            [ head ] @ tail
            |> List.map (fun { Line = (p1, p2); AtomCentersPosition = pos } ->
                let aDist = dist refTheta (direction p1)
                let bDist = dist refTheta (direction p2)
                if aDist < bDist then
                    aDist, bDist, { Line = (p1, p2); AtomCentersPosition = pos }
                else
                    bDist, aDist, { Line = (p2, p1); AtomCentersPosition = pos })
            |> List.sortBy (fun (a1, _, _) -> a1)

        let sections =
            [ for section in sections do
                  let a1, a2, clipping = section
                  let eclipsed =
                      sections
                      |> List.map (fun (b1, b2, _) -> b1 < a1 && b2 > a2)
                      |> List.contains true
                  match eclipsed with
                  | true -> ()
                  | false -> section ]

        let mutable s: Point2D option = None
        let rec reconstruct merged toMerge : Path =
            match toMerge with
            | [] -> merged
            | [c] ->
                let _, _, { Line = (p1, p2); AtomCentersPosition = pos } = c
                let largeArcFlag, sweepFlag = determineFlags p1 p2 pos
                match s with
                | None -> [ start p1; line p2; arc p1 largeArcFlag sweepFlag ]
//                | Some s -> merged @ [ arc p1 largeArcFlag sweepFlag; line p2; arc s largeArcFlag sweepFlag ]
                | Some s -> merged @ [ arc p1 Zero One; line p2; arc s Zero One ]
            | c1::c2::tail ->
                let a1a, a1b, { Line = (p1a, p1b); AtomCentersPosition = pos1 } = c1
                let largeArcFlag1, sweepFlag1 = determineFlags p1a p1b pos1
                let a2a, a2b, { Line = (p2a, p2b); AtomCentersPosition = pos2 } = c2
                let largeArcFlag2, sweepFlag2 = determineFlags p2a p2b pos2
                match a2a < a1b with
                | true ->
                    /// TODO: merge clippings; skip c2 for now
                    reconstruct merged ([ c1 ] @ tail)
                | false ->
                    match s with
                    | None ->
                        s <- Some p1a
                        reconstruct (merged @ [ start p1a; line p1b ]) ([ c2 ] @ tail)
                    | _ ->
//                        reconstruct (merged @ [ arc p1a largeArcFlag1 sweepFlag1; line p1b ]) ([ c2 ] @ tail)
                        reconstruct (merged @ [ arc p1a Zero One; line p1b ]) ([ c2 ] @ tail)

        [ opening ] @ (reconstruct [] sections) @ [ closing ]
        |> pathToSvg

let writeAtomsFilled (atoms: ProjectedAtomInfo list) : string =
    atoms
    |> List.map (fun atom -> drawAtom atom)
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
