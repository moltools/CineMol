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

let cinemolApi =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.buildProxy<ICinemolApi>

let init () : Model * Cmd<Msg> =
    let model =
        { Assignment = { Sdf = ""; Settings = { FilterHydrogens = false } }
          Svg = ""
          Encoded = "" }
    model, Cmd.none

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | UploadSdf (_, content) -> { model with Assignment = { model.Assignment with Sdf = content } }, Cmd.ofMsg Render
    | Render -> model, Cmd.OfAsync.perform cinemolApi.render model.Assignment GotEncoding
    | GotEncoding encoded -> { model with Encoded = encoded }, Cmd.none

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
                    reader.readAsText file ) ] ]

let private uploadFileButton dispatch =
    Html.div [
        prop.className "action-button"
        prop.children [
            Bulma.button.a [
                button.isOutlined
                button.isFullWidth
                prop.children [
                    Html.span "select SDF (Mol V2000) file"
                    uploadFileEvent dispatch ] ] ] ]

let private filterHydrogensSwitch dispatch =
    Bulma.field [
        Switch.checkbox [
            prop.id "mycheck"
            color.isDanger ]
        Html.label [
            prop.htmlFor "mycheck"
            prop.text "Check me" ] ]

let private svgViewer model dispatch =
    match model.Encoded with
    | e when e.Length = 0 -> Html.div [ prop.className "viewer" ]
    | _ ->
        Html.div [
            prop.className "viewer"
            prop.children [
                img [ Src $"data:image/svg+xml;base64,{model.Encoded}"] ] ]

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "cinemol"
        prop.children [
            uploadFileButton dispatch
            filterHydrogensSwitch dispatch
            svgViewer model dispatch ] ]