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

        // Parse atoms and apply rotation to point cloud
        let atoms =
            parseSdf sdf
            |> (fun atoms -> if showHydrogenAtoms = false then filterAtoms H atoms else atoms)
            |> Array.map (fun atom -> atom.Rotate Y ((rotX / 100.0) * 2.0 * Math.PI))
            |> Array.map (fun atom -> atom.Rotate Z ((rotY / 100.0) * 2.0 * Math.PI))
            |> Array.map (fun atom -> atom.Rotate X ((rotZ / 100.0) * 2.0 * Math.PI))

        let offsetViewBox =
            let minimumOffset = 2.0
            match atoms |> Array.map (fun a -> a.OriginalCenter.Distance origin) |> Array.max |> (*) 2.0 |> round 0 with
            | x when x < minimumOffset -> minimumOffset | x -> x

        let viewBox =
            match viewBox with
            | None -> -offsetViewBox, -offsetViewBox, offsetViewBox * 2.0, offsetViewBox * 2.0
            | Some x -> x

        let focalLength: float = offsetViewBox
        let pov: Point = { X = 1E-5; Y = 1E-5; Z = focalLength }
        let distPovOrigin: float = pov.Distance origin

        // Sort drawing order point cloud based on distance point to pov
        let atoms = atoms |> Array.sortBy (fun atom -> atom.OriginalCenter.Distance pov |> (*) -1.0)

        // Recalculate radius points based on distance point to pov
        let atoms = atoms |> Array.map (fun atom ->
            let projectedRadius = (distPovOrigin / (pov.Distance atom.OriginalCenter)) * atom.OriginalRadius
            { atom with ProjectedRadius = projectedRadius })

        let cameraForward: Vector = { X = -pov.X; Y = -pov.Y; Z = -pov.Z }
        let cameraPerpendicular: Vector = { X = cameraForward.Y; Y = -cameraForward.X; Z = 0.0 }
        let cameraHorizon: Vector = cameraForward.Cross cameraPerpendicular
        let setPerspective (atom: AtomInfo) : AtomInfo =
            let projectedCenter = project cameraPerpendicular cameraHorizon cameraForward (pov: Point) focalLength atom.OriginalCenter
            { atom with ProjectedCenter = projectedCenter }

        // Apply perspective projection on 3D point cloud on 2D view box
        let atoms = atoms |> Array.map (fun atom -> setPerspective atom)

        writeSVG viewBox depiction atoms, viewBox