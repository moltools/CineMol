module Index

open Browser
open Browser.Types
open Elmish
open Fable
open Fable.React
open Fable.React.Props
open Fable.Remoting.Client
open Feliz
open Feliz.Bulma
open Fulma
open Shared


type Model =
    { Assignment: Assignment
      Svg: string
      Encoded: string }

type Msg =
    | UploadSdf of name : string * content : string
    | Render
    | GotEncoding of string
    | SetHydrogenSwitch
    | SetRotation of float

let cinemolApi =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.buildProxy<ICinemolApi>

let init () : Model * Cmd<Msg> =
    let model =
        { Assignment = { Sdf = ""; Settings = { ShowHydrogenAtoms = false; Rotation = 0.5 } }
          Svg = ""
          Encoded = "" }
    model, Cmd.none

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | UploadSdf (_, content) -> { model with Assignment = { model.Assignment with Sdf = content } }, Cmd.ofMsg Render
    | Render -> model, Cmd.OfAsync.perform cinemolApi.render model.Assignment GotEncoding
    | GotEncoding encoded -> { model with Encoded = encoded }, Cmd.none
    | SetHydrogenSwitch ->
        let switch = not model.Assignment.Settings.ShowHydrogenAtoms
        let assignment = { model.Assignment with Settings = { model.Assignment.Settings with ShowHydrogenAtoms = switch } }
        { model with Assignment = assignment },
        Cmd.OfAsync.perform cinemolApi.render assignment GotEncoding
    | SetRotation f ->
        { model with Assignment = { model.Assignment with Settings = { model.Assignment.Settings with Rotation = f } } },
        Cmd.OfAsync.perform cinemolApi.render model.Assignment GotEncoding

let private uploadFileEvent dispatch =
    Fulma.File.input [
        Props [
            OnInput (
                fun ev ->
                    let file = (ev.target :?> HTMLInputElement).files.[0]
                    let name = file.name
                    let reader = FileReader.Create()
                    reader.onload <- fun _ ->
                        let content = reader.result :?> string
                        (name, content) |> UploadSdf |> dispatch
                    reader.readAsText file )
        ]
    ]

let private uploadFileButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                button.isOutlined
                button.isFullWidth
                color.isBlack
                prop.children [
                    Html.span "select SDF (Mol V2000) file"
                    uploadFileEvent dispatch
                ]
            ]
        ]
    ]

let private showHydrogenSwitch dispatch =
    Html.div [
        Switch.checkbox [
            prop.id "hydrogen-switch"
            prop.onChange (fun (_ : Event) -> SetHydrogenSwitch |> dispatch)
            color.isDanger ]
        Html.label [
            prop.htmlFor "hydrogen-switch"
            prop.text "Show hydrogen atoms"
        ]
    ]

let private svgViewer model =
    let svg =
        match model.Encoded with
        | e when e.Length = 0 -> ""
        | _ -> $"data:image/svg+xml;base64,{model.Encoded}"

    Html.div [
        prop.className "viewer"
        prop.children [
            img [
                Class "svg"
                Src svg
            ]
        ]
    ]

let private rotateMolSlider dispatch =
    Html.div [
        Slider.slider [
            slider.isFullWidth
            slider.isCircle
            slider.isLarge
            color.isBlack
            prop.onChange (fun (ev: Event) -> (SetRotation (float ev.Value) |> dispatch))
        ]
    ]

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "cinemol"
        prop.children [
            uploadFileButton dispatch
            showHydrogenSwitch dispatch
            svgViewer model
            rotateMolSlider dispatch
        ]
    ]