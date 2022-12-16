module Client.CineMol.Drawing

open System

open Helpers
open Types
open Svg

type DrawOptions = {
    Depiction: Depiction
    ShowHydrogenAtoms: bool
}
    with
    static member init = {
            Depiction = Filled
            ShowHydrogenAtoms = false
        }

let filterAtoms (atomType: Atom) (atoms: AtomInfo[]) : AtomInfo[] =
    Array.filter(fun (atom: AtomInfo) -> atom.Type <> atomType) atoms

let rotateAtoms (axis: Axis) (rads: float) (atoms: AtomInfo[]) : AtomInfo[] =
    atoms
    |> Array.map (fun (atom: AtomInfo) ->
        atom.Rotate axis ((rads / 100.0) * 2.0 * Math.PI))

let changeDistanceToAtoms (ratio: float) (atoms: AtomInfo[]) : AtomInfo[] =
    atoms
    |> Array.map (fun a -> { a with
                                OriginalCenter = {
                                    X = a.OriginalCenter.X * ratio
                                    Y = a.OriginalCenter.Y * ratio
                                    Z = a.OriginalCenter.Z * ratio
                                }
                                OriginalRadius = a.OriginalRadius * ratio
                            })

let draw
    (viewBox: ViewBox option)
    (options: DrawOptions)
    (rotation: Rotation)
    (zoom: Zoom)
    (mol: Molecule)
    : string * ViewBox =

    /// Rotate atoms based on supplied rotations in radians around axes.
    let mol = {
        mol with Atoms = mol.Atoms
                         |> rotateAtoms Y rotation.AxisX
                         |> rotateAtoms Z rotation.AxisY
                         |> rotateAtoms X rotation.AxisZ
    }

    /// Calculate view box offset (set before zoom, otherwise view box changes with zoom).
    let offsetViewBox =
        let minimumOffset = 2.0
        match mol.Atoms |> Array.map (fun a -> a.OriginalCenter.Distance origin) |> Array.max |> (*) 2.0 |> round 0 with
        | x when x < minimumOffset -> minimumOffset | x -> x

    let viewBox =
        match viewBox with
        | None -> -offsetViewBox, -offsetViewBox, offsetViewBox * 2.0, offsetViewBox * 2.0
        | Some x -> x

    let focalLength: float = offsetViewBox
    let pov: Point = { X = 1E-5; Y = 1E-5; Z = focalLength }
    let distPovOrigin: float = pov.Distance origin

    /// Filter atoms based on atom type.
    let mol = { mol with Atoms = if options.ShowHydrogenAtoms = false then filterAtoms H mol.Atoms else mol.Atoms }

    /// Apply zoom based on supplied scroll.
    /// Distance to (0, 0) and size of items will proportionally shrink or grow
    /// based on scroll direction.
    let mol = { mol with Atoms = changeDistanceToAtoms zoom.Ratio mol.Atoms }

    /// Sort drawing order point cloud based on distance point to POV.
    let mol = { mol with Atoms = mol.Atoms |> Array.sortBy (fun atom -> atom.OriginalCenter.Distance pov |> (*) -1.0) }

    /// Recalculate radius points based on distance point to POV.
    let mol = { mol with Atoms = mol.Atoms
                                 |> Array.map (fun atom ->
                                     let projectedRadius = (distPovOrigin / (pov.Distance atom.OriginalCenter)) * atom.OriginalRadius
                                     { atom with ProjectedRadius = projectedRadius })
    }

    let cameraForward: Vector = { X = -pov.X; Y = -pov.Y; Z = -pov.Z }
    let cameraPerpendicular: Vector = { X = cameraForward.Y; Y = -cameraForward.X; Z = 0.0 }
    let cameraHorizon: Vector = cameraForward.Cross cameraPerpendicular
    let setPerspective (atom: AtomInfo) : AtomInfo =
        let projectedCenter = project cameraPerpendicular cameraHorizon cameraForward (pov: Point) focalLength atom.OriginalCenter
        { atom with ProjectedCenter = projectedCenter }

    /// Apply perspective projection on 3D point cloud on 2D view box
    let mol = { mol with Atoms = mol.Atoms |> Array.map (fun atom -> setPerspective atom) }

    writeSVG viewBox options.Depiction mol.Atoms, viewBox
