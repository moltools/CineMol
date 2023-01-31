module Client.Index

open Browser
open Browser.Types
open CineMol.Types
open Elmish
open Fable.Core
open Fable.React
open Fable.React.Props
open Feliz
open Feliz.Bulma
open Fulma

open CineMol.Types.Style
open CineMol.Parsing
open CineMol.Drawing

let render () =
    async {
        return Some "", Some ""
    }

type MousePosition = { X: float; Y: float }
type WheelPosition = { Delta: float }
type DragTarget = | Dragging | NoTarget

[<RequireQualifiedAccess>]
module Cmd =
    let ups messageCtor =
        let handler dispatch =
            window.addEventListener("mouseup", fun _ ->
                dispatch messageCtor)
        [ handler ]

    let move messageCtor =
        let handler dispatch =
            window.addEventListener("mousemove", fun ev ->
                let ev = ev :?> MouseEvent
                { X = ev.pageX; Y = ev.pageY } |> messageCtor |> dispatch)
        [ handler ]

    let wheel messageCtor =
        let handler dispatch =
            window.addEventListener("mousewheel", fun ev ->
                let ev = ev :?> MouseWheelEvent
                { Delta = ev.wheelDelta } |> messageCtor |> dispatch)
        [ handler ]

type ViewerBackgroundColor = | Black | White
    with
    override this.ToString() =
        match this with | Black -> "#000000" | White -> "#FFFFFF"

type Model =
    { SdfString: string option
      SvgString: string option
      EncodedSvgString: string option
      DragTarget: DragTarget
      ViewerBackgroundColor: ViewerBackgroundColor
      Depiction: Depiction
      ShowHydrogenAtoms: bool
      RotationOverXAxis: float
      RotationOverYAxis: float
      SidebarCollapsed: bool }
    with
    static member Init () =
        { SdfString = None
          SvgString = None
          EncodedSvgString = None
          DragTarget = NoTarget
          ViewerBackgroundColor = White
          Depiction = BallAndStick
          ShowHydrogenAtoms = false
          RotationOverXAxis = 0.0
          RotationOverYAxis = 0.0
          SidebarCollapsed = false }
    member this.Reset () =
        { Model.Init() with
            ViewerBackgroundColor = this.ViewerBackgroundColor
            SidebarCollapsed = this.SidebarCollapsed }

type Msg =
    // User interface.
    | UploadSdf of name: string * content: string
    | ResetViewer
    | DownloadSvg
    | ToggleShowHydrogenAtoms
    | ToggleDepiction
    | ToggleBackgroundColor
    | ToggleSidebar

    // Render SVG.
    | Render
    | GotEncodedSvg of svgString: string option * encodedSvgString: string option
    | SetRotation of MousePosition

    // Mouse events.
    | MouseUp
    | MouseMove of MousePosition
    | MouseDrag of MousePosition
    | MouseDragStarted of MousePosition
    | MouseDragEnded
    | MouseWheelScroll of WheelPosition

let init () : Model * Cmd<Msg> =
    Model.Init(),
    Cmd.batch [ Cmd.ups MouseUp; Cmd.move MouseMove; Cmd.wheel MouseWheelScroll ]

let downloadSvgEvent (svgString: string option) (fileName: string) =
    match svgString with
    | Some svgString ->
        let anchor = Dom.document.createElement "a"
        let encodedSvgString =
            svgString
            |> sprintf "data:text/plain;charset=utf-8,%s.svg"
            |> JS.encodeURI
            |> (fun msg -> msg.Replace("#", "%23"))
        anchor.setAttribute("href", encodedSvgString)
        anchor.setAttribute("download", fileName)
        anchor.click()
    | None -> ()

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | UploadSdf (_, content) ->
        { model with SdfString = Some content }, Cmd.ofMsg Render

    | ResetViewer ->
        let newModel = model.Reset()
        newModel, Cmd.OfAsync.perform render () GotEncodedSvg

    | DownloadSvg ->
        downloadSvgEvent model.SvgString "model"
        model, Cmd.none

    | ToggleShowHydrogenAtoms ->
        let newModel = { model with ShowHydrogenAtoms = not model.ShowHydrogenAtoms }
        newModel, Cmd.OfAsync.perform render () GotEncodedSvg

    | ToggleDepiction ->
        let depiction =
            if model.Depiction = SpaceFilling then BallAndStick
            elif model.Depiction = BallAndStick then WireFrame
            else SpaceFilling
        let newModel = { model with Depiction = depiction }
        newModel, Cmd.OfAsync.perform render () GotEncodedSvg

    | ToggleBackgroundColor ->
        let backgroundColor = if model.ViewerBackgroundColor = Black then White else Black
        let newModel = { model with ViewerBackgroundColor = backgroundColor }
        newModel, Cmd.none

    | ToggleSidebar ->
        let newModel = { model with SidebarCollapsed = not model.SidebarCollapsed }
        newModel, Cmd.none

    | Render -> model, Cmd.OfAsync.perform render () GotEncodedSvg

    | GotEncodedSvg (svgString, encodedSvgString) ->
        { model with SvgString = svgString; EncodedSvgString = encodedSvgString }, Cmd.none

    | SetRotation position ->
        { model with RotationOverXAxis = position.X; RotationOverYAxis = position.Y },
        Cmd.OfAsync.perform render () GotEncodedSvg

    | MouseUp -> model, Cmd.ofMsg MouseDragEnded

    | MouseMove (position: MousePosition) -> model, Cmd.ofMsg (MouseDrag position)

    | MouseDragStarted _ -> { model with DragTarget = Dragging }, Cmd.none

    | MouseDragEnded -> { model with DragTarget = NoTarget }, Cmd.none

    | MouseDrag (position: MousePosition) ->
        let cmd = match model.DragTarget with | Dragging -> Cmd.ofMsg (SetRotation position) | _ -> Cmd.none
        model, cmd

    | MouseWheelScroll (position: WheelPosition) ->
        // TODO
        model, Cmd.none

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
                    reader.readAsText file

                    // Setting value to empty string removes cached filename again and makes sure you can upload the
                    // same file (after refresh) for a second time.
                    (ev.target :?> HTMLInputElement).value <- ""
            )
        ]
    ]

let sidebar model dispatch =
        let sidebarButton (action, icon, label: string) =
            Bulma.button.a [
                prop.className ("sidebar-button " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
                prop.children [
                    Html.i [ prop.className $"fas {icon}" ]
                    Html.span [ prop.style [ style.marginLeft (length.em 0.5) ]; prop.text (if model.SidebarCollapsed then "" else label) ]
                ]
                prop.onClick (fun _ -> dispatch action)
            ]

        Html.div [
            prop.className ("sidebar " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
            [ (ToggleBackgroundColor, "fa-adjust", "Toggle background")
              (ResetViewer, "fa-sync", "Refresh")
              (DownloadSvg, "fa-download", "Download")
              (ToggleShowHydrogenAtoms, "", "Toggle hydrogens")
              (ToggleDepiction, "fa-eye", "Toggle depiction") ]
            |> List.map (fun args -> sidebarButton args)
            |> prop.children
        ]

let viewer model dispatch =
    // let svg = draw()
    Html.div []

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "viewer"
        prop.style [ style.backgroundColor (model.ViewerBackgroundColor.ToString()) ]
        prop.children [
            sidebar model dispatch
            viewer model dispatch
        ]
    ]