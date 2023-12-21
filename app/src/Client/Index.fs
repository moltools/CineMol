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

open CineMol
open CineMol.Encoding
open CineMol.Types
open CineMol.Types.Geometry
open CineMol.Types.Chem
open CineMol.Types.Drawing
open CineMol.Parsing
open CineMol.Drawing

/// <summary>
/// Load example molecule penicillin G.
/// </summary>
let loadExampleMolecule () : string = """5904
  -OEChem-12182318453D

 41 43  0     1  0  0  0  0  0999 V2000
   -0.8019    1.2308    0.5170 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2842   -2.5451   -1.2026 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3517    1.0760   -0.8170 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1157   -0.6970    0.5961 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1598   -2.0405    1.2167 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4781   -0.7369    0.3018 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5677   -1.3807   -0.3375 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3100   -0.4177    1.1279 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.5670    1.6679    0.1134 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1671    0.3360   -0.3862 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.6142   -1.6325    0.4990 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.9181   -1.8261   -0.3007 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2193    2.2224    1.3871 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5771    2.7313   -0.9863 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6296    0.1576   -0.1297 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8730   -1.6107    0.1024 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9020   -1.2655   -0.9563 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8679   -0.1949   -0.5225 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5670    1.1373   -0.7687 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0457   -0.5556    0.1171 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4619    2.1290   -0.3670 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9405    0.4362    0.5190 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6486    1.7786    0.2769 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5167   -0.4852    2.1999 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9961    0.1926   -1.4616 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4509   -2.4585    1.2025 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2655    1.4761    2.1882 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6724    3.0921    1.7709 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2446    2.5572    1.1944 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0666    3.6467   -0.6644 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0842    2.3779   -1.8995 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6008    3.0119   -1.2561 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4359   -1.0113   -1.2762 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3150    0.9762   -0.6611 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4453   -2.1891   -1.1950 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4229   -0.9646   -1.8960 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6443    1.4206   -1.2675 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818   -1.5979    0.3120 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2343    3.1743   -0.5548 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8642    0.1634    1.0209 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.3451    2.5508    0.5901 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  8  1  0  0  0  0
  1  9  1  0  0  0  0
  2 12  2  0  0  0  0
  3 15  1  0  0  0  0
  3 34  1  0  0  0  0
  4 15  2  0  0  0  0
  5 16  2  0  0  0  0
  6  8  1  0  0  0  0
  6 10  1  0  0  0  0
  6 12  1  0  0  0  0
  7 11  1  0  0  0  0
  7 16  1  0  0  0  0
  7 33  1  0  0  0  0
  8 11  1  0  0  0  0
  8 24  1  0  0  0  0
  9 10  1  0  0  0  0
  9 13  1  0  0  0  0
  9 14  1  0  0  0  0
 10 15  1  0  0  0  0
 10 25  1  0  0  0  0
 11 12  1  0  0  0  0
 11 26  1  0  0  0  0
 13 27  1  0  0  0  0
 13 28  1  0  0  0  0
 13 29  1  0  0  0  0
 14 30  1  0  0  0  0
 14 31  1  0  0  0  0
 14 32  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 17 35  1  0  0  0  0
 17 36  1  0  0  0  0
 18 19  2  0  0  0  0
 18 20  1  0  0  0  0
 19 21  1  0  0  0  0
 19 37  1  0  0  0  0
 20 22  2  0  0  0  0
 20 38  1  0  0  0  0
 21 23  2  0  0  0  0
 21 39  1  0  0  0  0
 22 23  1  0  0  0  0
 22 40  1  0  0  0  0
 23 41  1  0  0  0  0
M  END
$$$$
"""

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
        match this with | Dark -> "#000" | Light -> "#f0f0f0"

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
      PreviousMousePosition: MousePosition option
      IsRendering: bool }
    with
    static member New () =
        { Molecule = None
          DrawingOptions = DrawingOptions.New()
          SvgString = None
          EncodedSvgString = None
          DragTarget = NoTarget
          ViewerBackgroundColor = Light
          SidebarCollapsed = false
          PreviousMousePosition = None
          IsRendering = false }
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
    | ToggleArt
    | ToggleBackgroundColor
    | ToggleHydrogenAtoms
    | LoadExample
    | HighResolution
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
        let depiction : ModelStyle =
            match model.DrawingOptions.ModelStyle with
            | SpaceFilling -> BallAndStick
            | BallAndStick -> WireFrame
            | WireFrame -> SpaceFilling
        let newModel = { model with DrawingOptions = { model.DrawingOptions with ModelStyle = depiction } }
        newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Toggle art style.
    | ToggleArt ->
        let art : ArtStyle =
            match model.DrawingOptions.ArtStyle with
            | Cartoon -> Glossy
            | Glossy -> Cartoon
        let newModel = { model with DrawingOptions = { model.DrawingOptions with ArtStyle = art } }
        newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Toggle background color (i.e., dark and light mode).
    | ToggleBackgroundColor ->
        let backgroundColor = if model.ViewerBackgroundColor = Dark then Light else Dark
        let newModel = { model with ViewerBackgroundColor = backgroundColor }
        newModel, Cmd.none

    // Toggle displaying of hydrogen atoms.
    | ToggleHydrogenAtoms ->
        let toggle = not model.DrawingOptions.DisplayHydrogenAtoms
        let newModel = { model with DrawingOptions = { model.DrawingOptions with DisplayHydrogenAtoms = toggle } }
        newModel, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Load  example molecule.
    | LoadExample ->
        let content = loadExampleMolecule()
        match (FileParser Sdf).Parse content with
        | Some (molecule::_) -> { model with Molecule = Some molecule }, Cmd.ofMsg Render
        | _ ->
            printfn "ERROR: Could not parse example molecule."
            model, Cmd.none // TODO: get error toast when parsing fails/file is empty

    // Draw current model in high resolution.
    | HighResolution ->
        let newModel = { model with DrawingOptions = { model.DrawingOptions with Resolution = 200 } }
        { newModel with IsRendering = true }, Cmd.OfAsync.perform render newModel GotEncodedSvg

    // Toggle sidebar between expanded and collapsed.
    | ToggleSidebar ->
        let newModel = { model with SidebarCollapsed = not model.SidebarCollapsed }
        newModel, Cmd.none

    // 'Render' molecule model.
    | Render -> model, Cmd.OfAsync.perform render model GotEncodedSvg

    // Process 'rendered' model.
    | GotEncodedSvg (svgString, encodedSvgString) ->
        let defaultOptions = DrawingOptions.New()
        { model with
            DrawingOptions = { model.DrawingOptions with Resolution = defaultOptions.Resolution }
            IsRendering = false
            SvgString = svgString
            EncodedSvgString = encodedSvgString }
        , Cmd.none

    // Set new rotation of molecule model.
    | SetRotation (x, y) ->
        match model.Molecule with
        | Some molecule ->
            let rotateAtom (atom : Atom) =
                let rotate axis pos (p: Point3D) = p.Rotate axis pos
                let rotatedCenter = atom.Position |> rotate Axis.X x |> rotate Axis.Y y
                { atom with Position = rotatedCenter }
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
                    Cmd.ofMsg (SetRotation (xRotation, -1.0 * yRotation))
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

    let rec toggleHydrogenAtoms dispatch =
        Bulma.button.a [
            prop.className ("sidebar-button " + match model.SidebarCollapsed with | true -> "collapsed" | false -> "expanded")
            prop.children [
                Fa.i [ Fa.Solid.Atom ] []
                Html.span [ prop.style [ style.marginLeft (length.em 0.5) ]; prop.text (if model.SidebarCollapsed then "" else "Toggle hydrogens") ]
            ]
            prop.onClick (fun _ -> dispatch ToggleHydrogenAtoms)
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

            Html.span [
                prop.className "version"
                if model.SidebarCollapsed then prop.children [ Html.div [ prop.className "version collapsed" ] ]
                else prop.text $"CineMol v{Version.version}"
            ]

            uploadFileButton dispatch
            generalButton DownloadSvg Fa.Solid.Download "Download"
            generalButton ResetViewer Fa.Solid.Sync "Refresh"
            generalButton ToggleBackgroundColor Fa.Solid.Adjust "Toggle background"
            generalButton ToggleDepiction Fa.Solid.Eye "Toggle depiction"
            generalButton ToggleArt Fa.Solid.PaintBrush "Toggle styling"
            toggleHydrogenAtoms dispatch
            generalButton LoadExample Fa.Solid.Seedling "Load example"
            generalButton HighResolution Fa.Solid.DrawPolygon "High resolution"
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
/// Loading indicator.
/// </summary>
let viewerRendering =
    Html.div [
        prop.className "loading"
        prop.children [
            Html.span [ prop.text "Loading..." ]
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
            match model.IsRendering with
            | false -> viewer model dispatch
            | true -> viewerRendering
        ]
    ]