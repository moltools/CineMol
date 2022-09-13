module Client.Index

open Browser
open Browser.Types
open Elmish
open Fable
open Fable.Core
open Fable.React
open Fable.React.Props
open Fable.Remoting.Client
open Feliz
open Feliz.Bulma
open Fulma
open Shared
open System


//type Position = { X: float; Y: float }

//type DragTarget =
//    | Dragging
//    | NoTarget

//[<RequireQualifiedAccess>]
//module Cmd =
//    let ups messageCtor =
//        let handler dispatch = window.addEventListener("mouseup", fun _ -> dispatch messageCtor)
//        [ handler ]
//
//    let move messageCtor =
//        let handler dispatch =
//            window.addEventListener("mousemove", fun ev ->
//                let ev = ev :?> MouseEvent
//                { X = ev.pageX; Y = ev.pageY } |> messageCtor |> dispatch)
//        [ handler ]

type Model =
    { Assignment: Assignment
      Svg: string
      Encoded: string
      Sidebar: Sidebar.Model }
//      DragTarget: DragTarget }

type Msg =
    | UploadSdf of name : string * content : string
    | Render
    | GotEncoding of svg : string * encodedSvg : string
    | SidebarMsg of Sidebar.Msg
    | SetRotation of float
    | MouseUp
//    | MouseMove of Position
//    | MouseDrag of Position
//    | MouseDragStarted of Guid * Position
    | MouseDragEnded

let cinemolApi =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.buildProxy<ICinemolApi>

let init () : Model * Cmd<Msg> =
    let sidebarModel, sidebarCmd = Sidebar.init()
    let cmd = Cmd.batch [
        Cmd.map SidebarMsg sidebarCmd
//        Cmd.ups MouseUp
//        Cmd.move MouseMove
    ]
    let model =
        { Assignment = { Sdf = ""; Settings = { ShowHydrogenAtoms = false; Rotation = 0.5 } }
          Svg = ""
          Encoded = ""
          Sidebar = sidebarModel }
//          DragTarget = NoTarget }
    model, cmd

let downloadSvg (svg : string) =
    let anchor = Dom.document.createElement "a"
    let contentReplace (oldValue: string) (newValue: string) (msg: string) = msg.Replace(oldValue, newValue)
    let encodedContent =
        svg
        |> sprintf "data:text/plain;charset=utf-8,%s"
        |> JS.encodeURI
        |> contentReplace "#" "%23"
    anchor.setAttribute("href", encodedContent)
    anchor.setAttribute("download", "model.svg")
    anchor.click()

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | UploadSdf (_, content) -> { model with Assignment = { model.Assignment with Sdf = content } }, Cmd.ofMsg Render
    | Render -> model, Cmd.OfAsync.perform cinemolApi.render model.Assignment GotEncoding
    | GotEncoding (svg, encodedSvg) -> { model with Encoded = encodedSvg; Svg = svg }, Cmd.none
    | SetRotation f ->
        { model with Assignment = { model.Assignment with Settings = { model.Assignment.Settings with Rotation = f } } },
        Cmd.OfAsync.perform cinemolApi.render model.Assignment GotEncoding
//    | MouseUp ->
//          model, Cmd.ofMsg MouseDragEnded
//    | MouseMove (position: Position) ->
//          model, Cmd.ofMsg (MouseDrag position)
//    | MouseDragStarted (guid, position) ->
//        { model with DragTarget = Dragging }, Cmd.none
//    | MouseDragEnded ->
//        { model with DragTarget = NoTarget }, Cmd.none
//    | MouseDrag (position: Position) ->
//        match model.DragTarget with
//        | Dragging -> model, Cmd.ofMsg (SetRotation position.X)
//        | _ -> model, Cmd.none
    | SidebarMsg msg ->
        let subModel, cmd, externalMsg = Sidebar.update msg model.Sidebar
        let newModel, extraCmd =
            match externalMsg with
            | Sidebar.NoOp ->
                model, Cmd.none
            | Sidebar.Reset ->
                { model with Encoded = "" }, Cmd.none
            | Sidebar.UploadSdf (_, src) ->
                { model with Assignment = { model.Assignment with Sdf = src } }, Cmd.ofMsg Render
            | Sidebar.DownloadSvg ->
                printf $"svg: {model.Svg}"
                downloadSvg model.Svg
                model, Cmd.none
            | Sidebar.GotShowHydrogenAtoms ->
                let switch = not model.Assignment.Settings.ShowHydrogenAtoms
                let assignment = { model.Assignment with Settings = { model.Assignment.Settings with ShowHydrogenAtoms = switch } }
                { model with Assignment = assignment },
                Cmd.OfAsync.perform cinemolApi.render assignment GotEncoding
            | Sidebar.GotRotation rotation ->
                model, Cmd.ofMsg (SetRotation rotation)
        { newModel with Sidebar = subModel }, Cmd.batch [ Cmd.map SidebarMsg cmd; extraCmd ]


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
            Sidebar.view model.Sidebar (SidebarMsg >> dispatch)
            svgViewer model
        ]
    ]