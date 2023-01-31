module Client.Index

open Browser
open Browser.Types
open Elmish
open Fable.Core
open Fable.FontAwesome
open Fable.React
open Fable.React.Props
open Feliz
open Feliz.Bulma
open Fulma

open CineMol.Encoding
open CineMol.Types.Chem
open CineMol.Types.Drawing
open CineMol.Parsing
open CineMol.Drawing

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
    { Molecule: Molecule option
      DrawingOptions: DrawingOptions option
      SvgString: string option
      EncodedSvgString: string option
      DragTarget: DragTarget
      ViewerBackgroundColor: ViewerBackgroundColor
      Depiction: ModelStyle
      ShowHydrogenAtoms: bool
      RotationOverXAxis: float
      RotationOverYAxis: float
      SidebarCollapsed: bool }
    with
    static member New () =
        { Molecule = None
          DrawingOptions = None
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
        { Model.New() with
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
    Model.New(),
    Cmd.batch [ Cmd.ups MouseUp; Cmd.move MouseMove; Cmd.wheel MouseWheelScroll ]

let uploadFileEvent dispatch =
    Fulma.File.input [
        Props [
            OnInput (
                fun ev ->
                    let file = (ev.target :?> HTMLInputElement).files[0]
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

let private downloadSvgEvent (svgString: string option) (fileName: string) =
    match svgString with
    | Some svgString ->
        let anchor = Dom.document.createElement "a"
        let encodedSvgString =
            svgString
            |> sprintf "data:text/plain;charset=utf-8,%s"
            |> JS.encodeURI
            |> (fun msg -> msg.Replace("#", "%23"))
        anchor.setAttribute("href", encodedSvgString)
        anchor.setAttribute("download", fileName)
        anchor.click()
    | None -> ()

let render (model: Model) =
    async {
        match model.Molecule with
        | Some molecule ->
            let svg = draw molecule
            let svgString = svg.ToString()
            let encodedSvgString = ISO_8859_1.Encode svgString |> sprintf "%A" // TODO: finsh encoding
            return Some svgString, Some encodedSvgString
        | None -> return None, None
    }

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | UploadSdf (_, content) ->
        match (FileParser Sdf).Parse content with
        | Some (molecule::_) -> { model with Molecule = Some molecule }, Cmd.ofMsg Render
        | _ -> model, Cmd.none // TODO: get error toast when parsing fails/file is empty

    | ResetViewer ->
        let newModel = model.Reset()
        newModel, Cmd.OfAsync.perform render model GotEncodedSvg

    | DownloadSvg ->
        downloadSvgEvent model.SvgString "model"
        model, Cmd.none

    | ToggleShowHydrogenAtoms ->
        let newModel = { model with ShowHydrogenAtoms = not model.ShowHydrogenAtoms }
        newModel, Cmd.OfAsync.perform render model GotEncodedSvg

    | ToggleDepiction ->
        let depiction =
            if model.Depiction = SpaceFilling then BallAndStick
            elif model.Depiction = BallAndStick then WireFrame
            else SpaceFilling
        let newModel = { model with Depiction = depiction }
        newModel, Cmd.OfAsync.perform render model GotEncodedSvg

    | ToggleBackgroundColor ->
        let backgroundColor = if model.ViewerBackgroundColor = Black then White else Black
        let newModel = { model with ViewerBackgroundColor = backgroundColor }
        newModel, Cmd.none

    | ToggleSidebar ->
        let newModel = { model with SidebarCollapsed = not model.SidebarCollapsed }
        newModel, Cmd.none

    | Render -> model, Cmd.OfAsync.perform render model GotEncodedSvg

    | GotEncodedSvg (svgString, encodedSvgString) ->
        { model with SvgString = svgString; EncodedSvgString = encodedSvgString }, Cmd.none

    | SetRotation position ->
        { model with RotationOverXAxis = position.X; RotationOverYAxis = position.Y },
        Cmd.OfAsync.perform render model GotEncodedSvg

    | MouseUp -> model, Cmd.ofMsg MouseDragEnded

    | MouseMove (position: MousePosition) -> model, Cmd.ofMsg (MouseDrag position)

    | MouseDragStarted _ -> { model with DragTarget = Dragging }, Cmd.none

    | MouseDragEnded -> { model with DragTarget = NoTarget }, Cmd.none

    | MouseDrag (position: MousePosition) ->
        let cmd = match model.DragTarget with | Dragging -> Cmd.ofMsg (SetRotation position) | _ -> Cmd.none
        model, cmd

    | MouseWheelScroll (_: WheelPosition) ->
        // TODO
        model, Cmd.none

let sidebar model dispatch =
        let uploadFileButton dispatch =
            Bulma.button.a [
                prop.className ("sidebar-button " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
                prop.children [
                    Fa.i [ Fa.Solid.Upload ] []
                    Html.span [ prop.style [ style.marginLeft (length.em 0.5) ]; prop.text (if model.SidebarCollapsed then "" else "Upload file") ]
                    uploadFileEvent dispatch
                ]
            ]

        let generalButton action  icon  (label: string) =
            Bulma.button.a [
                prop.className ("sidebar-button " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
                prop.children [
                    Fa.i [ icon ] []
                    Html.span [ prop.style [ style.marginLeft (length.em 0.5) ]; prop.text (if model.SidebarCollapsed then "" else label) ]
                ]
                prop.onClick (fun _ -> dispatch action)
            ]

        let reportBugButton =
            Bulma.button.a [
                prop.className ("sidebar-button " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
                prop.href "https://github.com/moltools/cinemol/issues"
                prop.rel "noreffer noopener"
                prop.target "_blank"
                prop.children [
                    Fa.i [ Fa.Solid.Bug ] []
                    Html.span [ prop.style [ style.marginLeft (length.em 0.5) ]; prop.text (if model.SidebarCollapsed then "" else "Report bug") ]
                ]
            ]

        Html.div [
            prop.className ("sidebar " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
            prop.children [
                uploadFileButton dispatch
                generalButton DownloadSvg Fa.Solid.Download "Download"
                generalButton ResetViewer Fa.Solid.Sync "Refresh"
                generalButton ToggleBackgroundColor Fa.Solid.Adjust "Toggle background"
                generalButton ToggleShowHydrogenAtoms Fa.Solid.ToggleOn "Toggle hydrogens"
                generalButton ToggleDepiction Fa.Solid.Eye "Toggle depiction"
                reportBugButton

                Html.div [
                    prop.className "collapse-button"
                    prop.children [ Fa.i [ (if model.SidebarCollapsed then Fa.Solid.AngleDoubleRight else Fa.Solid.AngleDoubleLeft) ] [] ]
                    prop.onClick (fun _ -> ToggleSidebar |> dispatch)
                ]
            ]
        ]

let viewer model dispatch =
    let viewerStyle =
        let sidebarWidth = match model.SidebarCollapsed with | true -> 60 | false -> 210
        let width = int window.innerWidth - sidebarWidth
        let height = int window.innerHeight
        if width < height then [ style.width width; style.height width ]
        else [ style.width height; style.height height ]

    let encodedSvgString = match model.EncodedSvgString with | None -> "" | Some s -> $"data:image/svg+xml;base64,{s}"

    Html.div [
        prop.className "viewer"
        prop.style viewerStyle
        prop.onMouseDown (fun ev ->
            ev.preventDefault()
            let coordsMouseDown = { X= ev.pageX; Y = ev.pageY }
            dispatch (MouseDragStarted coordsMouseDown))
        prop.children [
            img [
                Class "viewer"
                Style [ BackgroundColor (model.ViewerBackgroundColor.ToString()) ]
                Src encodedSvgString
            ]
        ]
    ]

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "app"
        prop.style [ style.backgroundColor (model.ViewerBackgroundColor.ToString()) ]
        prop.children [
            sidebar model dispatch
            viewer model dispatch
        ]
    ]