module Client.CineMol.Render

open Helpers
open Drawing

type Depiction =
    | Filled
    | BallAndStick

type ViewBox = float * float * float * float

type Settings =
    { ViewBox: ViewBox option
      Depiction: Depiction
      ShowHydrogenAtoms: bool
      XRotation: float
      YRotation: float
      ZRotation: float }

type Assignment = { Settings: Settings; Sdf: string }

let render = fun assignment -> async {
    let depiction: Types.Depiction =
        match assignment.Settings.Depiction with
        | Filled -> Types.Filled
        | BallAndStick -> Types.BallAndStick

    let svg, viewBox : string * ViewBox =
        assignment.Sdf
        |> draw assignment.Settings.ViewBox
                depiction
                assignment.Settings.ShowHydrogenAtoms
                assignment.Settings.XRotation
                assignment.Settings.YRotation
                assignment.Settings.ZRotation

    let encodedSvg : string = svg |> toBase64String

    return (svg, encodedSvg, viewBox) }