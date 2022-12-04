module Client.CineMol.Drawing

open System

open Helpers
open Types
open Svg

/// <summary>
///     Drawing options.
/// </summary>
type DrawOptions = {
    Depiction: Depiction
    ShowHydrogenAtoms: bool
}
    with
    static member init = {
            Depiction = Filled
            ShowHydrogenAtoms = false
        }

/// <summary>
///     Filter atoms based on atom type.
/// </summary>
/// <param name="atomType">
///     Atom type to filter atoms on.
/// </param>
/// <param name="atoms">
///     Atoms to filter on atom type.
/// </param>
/// <returns>
///     Filtered atoms.
/// </returns>
let filterAtoms (atomType: Atom) (atoms: AtomInfo[]) : AtomInfo[] =
    Array.filter(fun (atom: AtomInfo) -> atom.Type <> atomType) atoms

/// <summary>
///     Rotate atoms around axis.
/// </summary>
/// <param name="axis">
///     Axis to rotate atoms around.
/// </param>
/// <param name="rads">
///     Number of radians to rotate atoms around axis.
/// </param>
/// <param name="atoms">
///     Atoms to rotate around axis.
/// </param>
/// <returns>
///     Rotated atoms.
/// </returns>
let rotateAtoms (axis: Axis) (rads: float) (atoms: AtomInfo[]) : AtomInfo[] =
    atoms
    |> Array.map (fun (atom: AtomInfo) ->
        atom.Rotate axis ((rads / 100.0) * 2.0 * Math.PI))

/// <summary>
///     Create SVG representation of molecule.
/// </summary>
/// <param name="viewBox">
///     View box for SVG. Will calcualte view box bounds based on dimensions
///     of molecule if not supplied.
/// </param>
/// <param name="options">
///     Drawing options.
/// </param>
/// <param name="rotation">
///     Rotation options.
/// </param>
/// <param name="mol">
///     Molecule to draw.
/// </param>
/// <returns>
///     Tuple of SVG string and current bounds view box.
/// </returns>
let draw
    (viewBox: ViewBox option)
    (options: DrawOptions)
    (rotation: Rotation)
    (mol: Molecule)
    : string * ViewBox =

    /// Filter atoms based on atom type.
    let atoms =
        if options.ShowHydrogenAtoms = false then
            filterAtoms H mol.Atoms
        else
            mol.Atoms

    /// Rotate atoms based on supplied rotations in radians around axes.
    let atoms =
        atoms
        |> rotateAtoms Y rotation.AxisX
        |> rotateAtoms Z rotation.AxisY
        |> rotateAtoms X rotation.AxisZ

    /// Calculate view box offset.
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

    /// Sort drawing order point cloud based on distance point to POV.
    let atoms = atoms |> Array.sortBy (fun atom -> atom.OriginalCenter.Distance pov |> (*) -1.0)

    /// Recalculate radius points based on distance point to POV.
    let atoms = atoms |> Array.map (fun atom ->
        let projectedRadius = (distPovOrigin / (pov.Distance atom.OriginalCenter)) * atom.OriginalRadius
        { atom with ProjectedRadius = projectedRadius })

    let cameraForward: Vector = { X = -pov.X; Y = -pov.Y; Z = -pov.Z }
    let cameraPerpendicular: Vector = { X = cameraForward.Y; Y = -cameraForward.X; Z = 0.0 }
    let cameraHorizon: Vector = cameraForward.Cross cameraPerpendicular
    let setPerspective (atom: AtomInfo) : AtomInfo =
        let projectedCenter = project cameraPerpendicular cameraHorizon cameraForward (pov: Point) focalLength atom.OriginalCenter
        { atom with ProjectedCenter = projectedCenter }

    /// Apply perspective projection on 3D point cloud on 2D view box
    let atoms = atoms |> Array.map (fun atom -> setPerspective atom)

    writeSVG viewBox options.Depiction atoms, viewBox
