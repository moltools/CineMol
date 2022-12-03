module Client.Index

open Browser
open Browser.Types
open Elmish
open Fable.Core
open Fable.React
open Fable.React.Props
open Feliz
open Feliz.Bulma
open Fulma

open Client.CineMol.Types
open Client.CineMol.Parsing
open Client.CineMol.Encoding
open Client.CineMol.Drawing

type Settings =
    { ViewBox: ViewBox option
      Depiction: Depiction
      ShowHydrogenAtoms: bool
      XRotation: float
      YRotation: float
      ZRotation: float }

type Assignment = { Settings: Settings; Sdf: string }

let render = fun assignment -> async {
    let svg, viewBox : string * ViewBox =
        assignment.Sdf
        |> parseSdf
        |> (fun mols -> mols.[0])  // Select first molecule from SDF.
        |> draw assignment.Settings.ViewBox
                { Depiction = assignment.Settings.Depiction; ShowHydrogenAtoms = assignment.Settings.ShowHydrogenAtoms }
                { AxisX = assignment.Settings.XRotation; AxisY = assignment.Settings.YRotation; AxisZ = assignment.Settings.ZRotation }

    let encodedSvg : string = svg |> toBase64String

    return (svg, encodedSvg, viewBox) }

type Position = { X: float; Y: float }

type DragTarget =
    | Dragging
    | NoTarget

[<RequireQualifiedAccess>]
module Cmd =
    let ups messageCtor =
        let handler dispatch = window.addEventListener("mouseup", fun _ -> dispatch messageCtor)
        [ handler ]

    let move messageCtor =
        let handler dispatch =
            window.addEventListener("mousemove", fun ev ->
                let ev = ev :?> MouseEvent
                { X = ev.pageX; Y = ev.pageY } |> messageCtor |> dispatch)
        [ handler ]

type Model =
    { Assignment: Assignment
      Svg: string
      Encoded: string
      Sidebar: Sidebar.Model
      DragTarget: DragTarget }

type Msg =
    | UploadSdf of name: string * content: string
    | Render
    | GotEncoding of svg: string * encodedSvg: string * viewBox: ViewBox
    | SidebarMsg of Sidebar.Msg
    | SetRotation of Position
    | MouseUp
    | MouseMove of Position
    | MouseDrag of Position
    | MouseDragStarted of Position
    | MouseDragEnded

let init () : Model * Cmd<Msg> =
    let sidebarModel, sidebarCmd = Sidebar.init()
    let cmd = Cmd.batch [
        Cmd.map SidebarMsg sidebarCmd
        Cmd.ups MouseUp
        Cmd.move MouseMove
    ]
    let model =
        { Assignment = { Sdf = ""; Settings = { ViewBox = None
                                                Depiction = Filled
                                                ShowHydrogenAtoms = false
                                                XRotation = 0.5
                                                YRotation = 0.5
                                                ZRotation = 0.5 } }
          Svg = ""
          Encoded = ""
          Sidebar = sidebarModel
          DragTarget = NoTarget }
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
    | UploadSdf (_, content) ->
        { model with Assignment = { model.Assignment with Sdf = content } },
        Cmd.ofMsg Render
    | Render ->
        model,
        Cmd.OfAsync.perform render model.Assignment GotEncoding
    | GotEncoding (svg, encodedSvg, viewBox) ->
        let assignment = { model.Assignment with Settings = { model.Assignment.Settings with ViewBox = Some viewBox } }
        { model with Encoded = encodedSvg; Svg = svg; Assignment = assignment },
        Cmd.none
    | SetRotation position ->
        { model with Assignment = { model.Assignment with Settings = { model.Assignment.Settings with XRotation = position.X; YRotation = position.Y } } },
        Cmd.OfAsync.perform render model.Assignment GotEncoding
    | MouseUp ->
          model, Cmd.ofMsg MouseDragEnded
    | MouseMove (position: Position) ->
          model, Cmd.ofMsg (MouseDrag position)
    | MouseDragStarted position ->
        { model with DragTarget = Dragging }, Cmd.none
    | MouseDragEnded ->
        { model with DragTarget = NoTarget }, Cmd.none
    | MouseDrag (position: Position) ->
        match model.DragTarget with
        | Dragging -> model, Cmd.ofMsg (SetRotation position)
        | _ -> model, Cmd.none
    | SidebarMsg msg ->
        let subModel, cmd, externalMsg = Sidebar.update msg model.Sidebar
        let newModel, extraCmd =
            match externalMsg with
            | Sidebar.NoOp ->
                model, Cmd.none
            | Sidebar.Reset ->
                let assignment = { model.Assignment with Sdf = ""; Settings = { model.Assignment.Settings with ViewBox = None } }
                { model with Encoded = ""; Svg = ""; Assignment = assignment },
                Cmd.OfAsync.perform render assignment GotEncoding
            | Sidebar.UploadSdf (_, src) ->
                let assignment = { model.Assignment with Sdf = src }
                { model with Assignment = assignment }, Cmd.ofMsg Render
            | Sidebar.DownloadSvg ->
                printf $"svg: {model.Svg}"
                downloadSvg model.Svg
                model, Cmd.none
            | Sidebar.GotShowHydrogenAtoms ->
                let switch = not model.Assignment.Settings.ShowHydrogenAtoms
                let assignment = { model.Assignment with Settings = { model.Assignment.Settings with ShowHydrogenAtoms = switch } }
                { model with Assignment = assignment },
                Cmd.OfAsync.perform render assignment GotEncoding
            | Sidebar.GotDepiction ->
                let switch = if model.Assignment.Settings.Depiction = BallAndStick then Filled else BallAndStick
                let assignment = { model.Assignment with Settings = { model.Assignment.Settings with Depiction = switch } }
                { model with Assignment = assignment },
                Cmd.OfAsync.perform render assignment GotEncoding
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

let private svgViewer (dispatch: Msg -> unit) model =
    let svg =
        match model.Encoded with
        | e when e.Length = 0 -> ""
        | _ -> $"data:image/svg+xml;base64,{model.Encoded}"

    Html.div [
        prop.className "viewer"
        prop.onMouseDown (fun ev ->
            ev.preventDefault()
            let coordsMouseDown = { X = ev.pageX; Y = ev.pageY }
            dispatch (MouseDragStarted coordsMouseDown))
        prop.children [
            img [
                Class "svg"
                Src svg
            ]
        ]
    ]

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "cinemole"
        prop.children [
            Sidebar.view model.Sidebar (SidebarMsg >> dispatch)
            svgViewer dispatch model
        ]
    ]