namespace Cinemol

module Drawing =

    open System

    open Cinemol.Helpers
    open Cinemol.Types
    open Cinemol.Parsing
    open Cinemol.Svg

    let draw
        (viewBox: ViewBox option)   // Needs to be set if None based on distances points in point cloud
        (depiction: Depiction)      // Filled or ball-and-stick
        (showHydrogenAtoms: bool)   // Filter out hydrogen atoms or not
        (rotX: float)               // Rotation around x-axis
        (rotY: float)               // Rotation around y-axis
        (rotZ: float)               // Rotation around z-axis
        (sdf: string)               // Input SDF (V2000) mol file in string format
        : string * ViewBox =        // Returns SVG string and (set) view box
        let atoms =
            parseSdf sdf
            |> (fun atoms -> if showHydrogenAtoms = false then filterAtoms H atoms else atoms)
            |> Array.map (fun atom -> atom.Rotate Y ((rotX / 100.0) * 2.0 * Math.PI))
            |> Array.map (fun atom -> atom.Rotate Z ((rotY / 100.0) * 2.0 * Math.PI))
            |> Array.map (fun atom -> atom.Rotate X ((rotZ / 100.0) * 2.0 * Math.PI))

        let offsetViewBox =
            let minimumOffset = 10.0
            match atoms |> Array.map (fun a -> a.Distance origin) |> Array.max |> (*) 2.0 |> round 0 with
            | x when x < minimumOffset -> minimumOffset | x -> x

        let viewBox =
            match viewBox with
            | None -> -offsetViewBox, -offsetViewBox, offsetViewBox * 2.0, offsetViewBox * 2.0
            | Some x -> x

        let pov: Point = { X = 0.0; Y = 0.0; Z = offsetViewBox * 2.0 }
        let cameraForward: Vector = { X = -pov.X; Y = -pov.Y; Z = -pov.Z }
        let cameraPerpendicular: Vector = { X = cameraForward.Y; Y = -cameraForward.X; Z = 0.0 }
        let cameraHorizon: Vector = cameraForward.Cross cameraPerpendicular
        let focalLength = 1000.0

        let warpAtom (atom: AtomInfo) : AtomInfo =
            { atom with Center = warp cameraPerpendicular cameraHorizon cameraForward (pov: Point) focalLength atom.Center }

        let atoms =
            atoms
            |> Array.map (fun atom -> warpAtom atom)
            |> Array.sortBy (fun atom -> atom.Center.Z)

        // TODO: fix bug -- all dot products in project vector are 0.0 and result in NaN when dividing by magnitude
        // TODO: apply perspective to atom coordinates and atom radius based on camera position
        // TODO: clipping

        writeSVG viewBox depiction atoms, viewBox