module Client.CineMol.Drawing

open System

open Helpers
open Styles
open Types
open Geometry
open Svg

type DrawOptions = {Depiction: Depiction; ShowHydrogenAtoms: bool }
    with static member init = { Depiction = Filled; ShowHydrogenAtoms = false }

let filterAtoms (atomType: AtomType) (atoms: AtomInfo[]) : AtomInfo[] =
    Array.filter(fun (atom: AtomInfo) -> atom.AtomType <> atomType) atoms

let rotateAtoms (axis: Axis) (rads: float) (atoms: AtomInfo[]) : AtomInfo[] =
    atoms
    |> Array.map (fun (atom: AtomInfo) ->
        atom.Rotate axis ((rads / 100.0) * 2.0 * Math.PI))

let changeDistanceToAtoms (ratio: float) (mol: Molecule) : Molecule =
    let transformAtom (a: AtomInfo) =
        { a with
            Center = {
                X = a.Center.X * ratio
                Y = a.Center.Y * ratio
                Z = a.Center.Z * ratio  }
            Radius = a.Radius * ratio }
    let transformBond b =
        { b with Scaling = b.Scaling * ratio }
    { mol with Atoms = Array.map (fun a -> transformAtom a) mol.Atoms; Bonds = Array.map (fun b -> transformBond b) mol.Bonds }

let draw
    (viewBox: ViewBox option)
    (options: DrawOptions)
    (rotation: Rotation)
    (zoom: Zoom)
    (mol: Molecule)
    : string * ViewBox =

    // Rotate atoms based on supplied rotations in radians around axes.
    let mol = {
        mol with Atoms = mol.Atoms
                         |> rotateAtoms Y rotation.AxisX
                         |> rotateAtoms Z rotation.AxisY
                         |> rotateAtoms X rotation.AxisY
    }

    // Calculate view box offset (set before zoom, otherwise view box changes with zoom).
    let offsetViewBox =
        let minimumOffset = 2.0
        match mol.Atoms |> Array.map (fun a -> a.Center.Distance origin) |> Array.max |> (*) 2.0 |> round 0 with
        | x when x < minimumOffset -> minimumOffset | x -> x

    let viewBox =
        match viewBox with
        | None -> -offsetViewBox, -offsetViewBox, offsetViewBox * 2.0, offsetViewBox * 2.0
        | Some x -> x

    let focalLength: float = offsetViewBox
    let pov: Point3D = { X = 1E-5; Y = 1E-5; Z = focalLength }
    let distPovOrigin: float = pov.Distance origin

    // Filter atoms based on atom type.
    let mol = { mol with Atoms = if options.ShowHydrogenAtoms = false then filterAtoms H mol.Atoms else mol.Atoms }

    // Apply zoom based on supplied scroll.
    let perspectiveMol = changeDistanceToAtoms zoom.Ratio mol

    // Sort drawing order point cloud based on distance point to POV.
    let perspectiveMol = { perspectiveMol with Atoms = perspectiveMol.Atoms |> Array.sortBy (fun atom -> atom.Center.Distance pov |> (*) -1.0) }

    // Recalculate radii based on distance point to POV.
    let perspectiveMol = {
        perspectiveMol with
            Atoms = perspectiveMol.Atoms
                    |> Array.map (fun atom ->
                        let projectedRadius = (distPovOrigin / (pov.Distance atom.Center)) * atom.Radius
                        { atom with Radius = projectedRadius }) }

    // Apply perspective projection on 3D point cloud on 2D view box
    let cameraForward: Vector3D = { X = -pov.X; Y = -pov.Y; Z = -pov.Z }
    let cameraPerpendicular: Vector3D = { X = cameraForward.Y; Y = -cameraForward.X; Z = 0.0 }
    let cameraHorizon: Vector3D = cameraForward.Cross cameraPerpendicular
    let project p = project cameraPerpendicular cameraHorizon cameraForward pov focalLength p
    let setPerspectiveAtom (atom: AtomInfo) : AtomInfo = { atom with Center = project atom.Center }

    let perspectiveMol = {
        perspectiveMol with
            Atoms = perspectiveMol.Atoms |> Array.map (fun atom -> setPerspectiveAtom atom)
    }

    // Calculate clipping.
    let projectedMol =
        { Atoms =
            Array.zip perspectiveMol.Atoms mol.Atoms
            |> Array.map (fun (perspectiveAtom, atom) ->
                { Index = perspectiveAtom.Index
                  AtomType = perspectiveAtom.AtomType
                  Center = { X = perspectiveAtom.Center.X; Y = perspectiveAtom.Center.Y }
                  Radius = perspectiveAtom.Radius
//                  Clipping = clip perspectiveAtom perspectiveMol atom mol
                  Clipping = []
                })
          Bonds = perspectiveMol.Bonds }

    writeSVG viewBox options.Depiction projectedMol, viewBox
