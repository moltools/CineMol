module Client.Index

open System

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
open CineMol.Types.Geometry
open CineMol.Types.Chem
open CineMol.Types.Drawing
open CineMol.Parsing
open CineMol.Drawing

/// <summary>
/// Position of mouse on screen.
/// </summary>
type MousePosition = { X: float; Y: float }

/// <summary>
/// Wheel scroll delta.
/// </summary>
type WheelPosition = { Delta: float }

/// <summary>
/// Dragging state.
/// </summary>
type DragTarget = | Dragging | NoTarget

/// <summary>
/// Create handlers for mouse events.
/// </summary>
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

/// <summary>
/// Viewer background colors (i.e., dark and light modes).
/// </summary>
type ViewerBackgroundColor = | Dark | Light
    with
    override this.ToString() =
        match this with | Dark -> "#343231" | Light -> "#F0F0F0"

/// <summary>
/// App model.
/// </summary>
type Model =
    { Molecule: Molecule option
      DrawingOptions: DrawingOptions
      SvgString: string option
      EncodedSvgString: string option
      DragTarget: DragTarget
      ViewerBackgroundColor: ViewerBackgroundColor
      SidebarCollapsed: bool
      PreviousMousePosition: MousePosition option }
    with
    static member New () =
        { Molecule = None
          DrawingOptions = { DrawingOptions.New() with Style = BallAndStick }
          SvgString = None
          EncodedSvgString = None
          DragTarget = NoTarget
          ViewerBackgroundColor = Dark
          SidebarCollapsed = false
          PreviousMousePosition = None }
    member this.Reset () =
        { Model.New() with
            ViewerBackgroundColor = this.ViewerBackgroundColor
            SidebarCollapsed = this.SidebarCollapsed }

/// <summary>
/// App messages.
/// </summary>
type Msg =
    // User interface.
    | UploadSdf of name: string * content: string
    | ResetViewer
    | DownloadSvg
    | ToggleDepiction
    | ToggleBackgroundColor
    | ToggleSidebar

    // Rendering.
    | Render
    | GotEncodedSvg of svgString: string option * encodedSvgString: string option
    | SetRotation of xRotation: float * yRotation: float

    // Mouse events.
    | MouseUp
    | MouseMove of MousePosition
    | MouseDrag of MousePosition
    | MouseDragStarted of MousePosition
    | MouseDragEnded
    | MouseWheelScroll of WheelPosition

/// <summary>
/// Initialize app.
/// </summary>
let init () : Model * Cmd<Msg> =
    Model.New(),
    Cmd.batch [ Cmd.ups MouseUp; Cmd.move MouseMove; Cmd.wheel MouseWheelScroll ]

/// <summary>
/// Element that creates file picker pop-up.
/// </summary>
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

/// <summary>
/// Download text to file to downloads folder.
/// </summary>
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

/// <summary>
/// 'Render' molecule model.
/// </summary>
let render (model: Model) =
    async {
        match model.Molecule with
        | Some molecule ->
            let svg, options = draw molecule model.DrawingOptions
            let svgString = svg.ToString()
            let encodedSvgString = ISO_8859_1.Encode svgString |> Convert.ToBase64String
            return Some svgString, Some encodedSvgString
        | None -> return None, None
    }

/// <summary>
/// Update model.
/// </summary>
let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with

    // Upload SDF file and read and parse contents into Molecule object.
    | UploadSdf (_, content) ->
        match (FileParser Sdf).Parse content with
        | Some (molecule::_) -> { model with Molecule = Some molecule }, Cmd.ofMsg Render
        | _ -> model, Cmd.none // TODO: get error toast when parsing fails/file is empty

    // Reset viewer.
    | ResetViewer ->
        let newModel = model.Reset()
        newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Download SVG to downloads folder.
    | DownloadSvg ->
        downloadSvgEvent model.SvgString "model.svg"
        model, Cmd.none

    // Toggle model depiction.
    | ToggleDepiction ->
        let depiction =
            if model.DrawingOptions.Style = SpaceFilling then BallAndStick
            elif model.DrawingOptions.Style = BallAndStick then Tube
            else SpaceFilling
        let newModel = { model with DrawingOptions = { model.DrawingOptions with Style = depiction } }
        newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Toggle background color (i.e., dark and light mode).
    | ToggleBackgroundColor ->
        let backgroundColor = if model.ViewerBackgroundColor = Dark then Light else Dark
        let newModel = { model with ViewerBackgroundColor = backgroundColor }
        newModel, Cmd.none

    // Toggle sidebar between expanded and collapsed.
    | ToggleSidebar ->
        let newModel = { model with SidebarCollapsed = not model.SidebarCollapsed }
        newModel, Cmd.none

    // 'Render' molecule model.
    | Render -> model, Cmd.OfAsync.perform render model GotEncodedSvg

    // Process 'rendered' model.
    | GotEncodedSvg (svgString, encodedSvgString) ->
        { model with SvgString = svgString; EncodedSvgString = encodedSvgString }, Cmd.none

    // Set new rotation of molecule model.
    | SetRotation (x, y) ->
        match model.Molecule with
        | Some molecule ->
            let rotateAtom (Atom3D(i, c, r)) =
                let rotate axis pos (p: Point3D) = p.Rotate axis pos
                let rotatedCenter = c |> rotate Axis.X x |> rotate Axis.Y y
                Atom3D (i, rotatedCenter, r)
            let rotatedAtoms = molecule.Atoms |> List.map (fun a -> rotateAtom a)
            let rotatedMolecule = { molecule with Atoms = rotatedAtoms }
            let newModel = { model with Molecule = Some rotatedMolecule }
            newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg
        | None -> model, Cmd.none

    // Mouse up indicates the end of a drag event.
    | MouseUp -> model, Cmd.ofMsg MouseDragEnded

    // Relay mouse drag.
    | MouseMove pos -> model, Cmd.ofMsg (MouseDrag pos)

    // Set start of drag event.
    | MouseDragStarted pos -> { model with DragTarget = Dragging }, Cmd.none

    // Set end of drag event.
    | MouseDragEnded -> { model with DragTarget = NoTarget }, Cmd.none

    // Register dragging movement.
    | MouseDrag (pos: MousePosition) ->
        let cmd =
            match model.PreviousMousePosition with
            | Some prevPos ->
                match model.DragTarget with
                | Dragging ->
                    let xRotation = (prevPos.Y - pos.Y) / 180.0
                    let yRotation = (prevPos.X - pos.X) / 180.0
                    Cmd.ofMsg (SetRotation (xRotation, yRotation))
                | _ -> Cmd.none
            | None -> Cmd.none

        let newModel = { model with PreviousMousePosition = Some pos }
        newModel, cmd

    // Register mouse scroll.
    | MouseWheelScroll (pos: WheelPosition) ->
        // TODO: implement zoom
        model, Cmd.none

/// <summary>
/// Sidebar element.
/// </summary>
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
                generalButton ToggleDepiction Fa.Solid.Eye "Toggle depiction"
                reportBugButton

                Html.div [
                    prop.className "collapse-button"
                    prop.children [ Fa.i [ (if model.SidebarCollapsed then Fa.Solid.AngleDoubleRight else Fa.Solid.AngleDoubleLeft) ] [] ]
                    prop.onClick (fun _ -> ToggleSidebar |> dispatch)
                ]
            ]
        ]

/// <summary>
/// Viewer element.
/// </summary>
let viewer model dispatch =
    let viewerStyle =
        let sidebarWidth = match model.SidebarCollapsed with | true -> 60 | false -> 210
        let width = int window.innerWidth - sidebarWidth
        let height = int window.innerHeight
        if width < height then [ style.width width; style.height width ]
        else [ style.width height; style.height height ]

    let encodedSvgString =
        match model.EncodedSvgString with
        | None -> ""
        | Some s -> sprintf "data:image/svg+xml;base64,%s" s

    Html.div [
        prop.className "viewer"
        prop.style viewerStyle
        prop.onMouseDown (fun ev ->
            ev.preventDefault()
            dispatch (MouseDragStarted { X = ev.pageX; Y = ev.pageY }))
        prop.children [
            img [ Style [ BackgroundColor (model.ViewerBackgroundColor.ToString()) ]; Src encodedSvgString ]
        ]
    ]

/// <summary>
/// Main app view.
/// </summary>
let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "app"
        prop.style [ style.backgroundColor (model.ViewerBackgroundColor.ToString()) ]
        prop.children [
            sidebar model dispatch
            viewer model dispatch
        ]
    ]