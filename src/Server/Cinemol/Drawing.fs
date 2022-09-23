namespace Cinemol

module Drawing =

    open System

    open Cinemol.Types
    open Cinemol.Parsing
    open Cinemol.Geometry
    open Cinemol.Svg

    let draw (viewBox: ViewBox option) depiction showHydrogenAtoms rotX rotY rotZ sdf =
        let filterAtoms (atomType : Atom) (atoms : AtomInfo array) : AtomInfo array =
            Array.filter(fun ((_, a, _) : AtomInfo) -> a <> atomType) atoms

        let atoms =
            match showHydrogenAtoms with
            | true -> parse_sdf sdf
            | false -> sdf |> parse_sdf |> filterAtoms H

        let maxX = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.X) atoms |> Array.max
        let maxY = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.Y) atoms |> Array.max
        let maxZ = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.Z) atoms |> Array.max
        let unitSize = (List.max [ maxX; maxY; maxZ ]) * 2.0
        let pov = { X = 0.0; Y = 0.0; Z = unitSize }

        let viewBox =
            match viewBox with
            | None ->
                let offset = unitSize |> int |> float
                - offset, - offset, offset * 2.0, offset * 2.0
            | Some vb -> vb

        let radX : float option = Some ((rotX / 100.0) * 2.0 * Math.PI)
        let radY : float option = Some ((rotY / 100.0) * 2.0 * Math.PI)
        let radZ : float option = Some ((rotZ / 100.0) * 2.0 * Math.PI)
        let rotatedAtoms =
            atoms
            |> Array.map (fun ((i, a, c) : AtomInfo) -> (i, a, rotateAxisY c radX))
            |> Array.map (fun ((i, a, c) : AtomInfo) -> (i, a, rotateAxisZ c radY))
            |> Array.map (fun ((i, a, c) : AtomInfo) -> (i, a, rotateAxisX c radZ))
            |> Array.sortBy (fun ((_, _, c) : AtomInfo) -> - (abs (pov.Z - c.Z)))

        writeSVG viewBox depiction pov unitSize rotatedAtoms, viewBox