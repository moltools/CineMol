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

// ============================================================================
// Rendering SVG.
// ============================================================================
let render =
    fun (vb: ViewBox option, options: DrawOptions, rot: Rotation, zoom: Zoom, sdf: string) ->
        async {
            let molToDraw =
                sdf
                |> parseSdf
                |> (fun mols -> mols.[0])

            let svg, vb: string * ViewBox option =
                if molToDraw.Atoms.Length = 0 then
                    "", None
                else
                    let svg, vb = draw vb options rot zoom molToDraw
                    svg, Some vb

            return (svg, svg |> toBase64String, vb)
        }

// ============================================================================
// Mouse events for SVG.
// ============================================================================
type MousePosition = { X: float; Y: float }

type WheelPosition = { Delta: float }

type DragTarget =
    | Dragging
    | NoTarget

[<RequireQualifiedAccess>]
module Cmd =
    let ups messageCtor =
        let handler dispatch =
            window.addEventListener("mouseup", fun _ -> dispatch messageCtor)
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

// ============================================================================
// Model.
// ============================================================================
type Model = {
    Sdf: string
    Svg: string
    EncodedSvg: string
    ViewBox: ViewBox option
    DrawOptions: DrawOptions
    Rotation: Rotation
    Zoom: Zoom
    DragTarget: DragTarget
}
    with
    member x.renderArgs = x.ViewBox, x.DrawOptions, x.Rotation, x.Zoom, x.Sdf
    static member init = {
            Sdf = ""
            Svg = ""
            EncodedSvg = ""
            DrawOptions = DrawOptions.init
            Rotation = Rotation.init
            Zoom = Zoom.init
            ViewBox = None
            DragTarget = NoTarget
        }

type Msg =
    /// GUI buttons.
    | UploadSdf of name: string * content: string
    | ResetViewer
    | DownloadSvg
    | ToggleShowHydrogenAtoms
    | ToggleDepiction

    /// Rendering SVG.
    | Render
    | GotEncoding of svg: string * encodedSvg: string * viewBox: ViewBox option
    | SetRotation of MousePosition

    /// Mouse events for dragging SVG.
    | MouseUp
    | MouseMove of MousePosition
    | MouseDrag of MousePosition
    | MouseDragStarted of MousePosition
    | MouseDragEnded

    /// Mouse events for zooming in and out on SVG.
    | WheelScroll of WheelPosition

// ============================================================================
// Actions.
// ============================================================================

/// <summary>
///     Custom action for downloading SVG.
/// </summary>
let downloadSvgEvent (svg : string) =
    let anchor = Dom.document.createElement "a"
    let contentReplace (oldValue: string) (newValue: string) (msg: string) =
        msg.Replace(oldValue, newValue)
    let encodedContent =
        svg
        |> sprintf "data:text/plain;charset=utf-8,%s"
        |> JS.encodeURI
        |> contentReplace "#" "%23"
    anchor.setAttribute("href", encodedContent)
    anchor.setAttribute("download", "model.svg")
    anchor.click()

// ============================================================================
// Initializing model.
// ============================================================================
let init () : Model * Cmd<Msg> =
    let cmd =
        Cmd.batch [
            Cmd.ups MouseUp
            Cmd.move MouseMove
            Cmd.wheel WheelScroll
        ]
    Model.init, cmd

// ============================================================================
// Updating model.
// ============================================================================
let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with

    /// Upload SDF.
    | UploadSdf (_, content) ->
        { model with Sdf = content }, Cmd.ofMsg Render

    /// Reset viewer.
    | ResetViewer ->
        let newModel = Model.init
        newModel, Cmd.OfAsync.perform render newModel.renderArgs GotEncoding

    /// Download SVG.
    | DownloadSvg ->
        downloadSvgEvent model.Svg
        model, Cmd.none

    /// Toggle show hydrogens.
    | ToggleShowHydrogenAtoms ->
        let toggle = not model.DrawOptions.ShowHydrogenAtoms

        let newModel = {
            model with
                DrawOptions = {
                    model.DrawOptions with
                        ShowHydrogenAtoms = toggle
                }
            }

        newModel, Cmd.OfAsync.perform render newModel.renderArgs GotEncoding

    /// Toggle depiction.
    | ToggleDepiction ->
        let toggle =
            if model.DrawOptions.Depiction = Depiction.Filled then
                Depiction.BallAndStick
            else
                Depiction.Filled

        let newModel = {
            model with
                DrawOptions = {
                    model.DrawOptions with
                        Depiction = toggle
                }
            }

        newModel, Cmd.OfAsync.perform render newModel.renderArgs GotEncoding

    /// Render SVG from SDF message.
    | Render ->
        model, Cmd.OfAsync.perform render model.renderArgs GotEncoding

    /// Store rendered, encoded SVG in model message.
    | GotEncoding (svg, encodedSvg, viewBox) ->
        { model with
            EncodedSvg = encodedSvg
            Svg = svg
            ViewBox = viewBox
        }, Cmd.none

    /// Set new rotation of molecule in SVG message.
    | SetRotation position ->
        { model with
            Rotation = {
                model.Rotation with
                    AxisX = position.X
                    AxisY = position.Y
            }
        },
        Cmd.OfAsync.perform render model.renderArgs GotEncoding

    /// Mouse drag messages for rendered, encoded SVG in model.
    | MouseUp ->
        model, Cmd.ofMsg MouseDragEnded

    | MouseMove (position: MousePosition) ->
        model, Cmd.ofMsg (MouseDrag position)

    | MouseDragStarted position ->
        { model with DragTarget = Dragging }, Cmd.none

    | MouseDragEnded ->
        { model with DragTarget = NoTarget }, Cmd.none

    | MouseDrag (position: MousePosition) ->
        match model.DragTarget with
        | Dragging -> model, Cmd.ofMsg (SetRotation position)
        | _ -> model, Cmd.none

    /// Mouse scroll message for rendered, encoded SVG in model.
    | WheelScroll (wheel: WheelPosition) ->
        let incr =
            if wheel.Delta < 0.0 then
                -0.05
            else
                +0.05

        let newRatio =
            let newRatio = model.Zoom.Ratio + incr
            if newRatio < 0.0 then 0.0 else newRatio 

        let newModel =
            { model with
                Zoom = {
                    model.Zoom with Ratio = newRatio
                }
            }
        newModel, Cmd.OfAsync.perform render newModel.renderArgs GotEncoding

// ============================================================================
// GUI element: upload file button.
// ============================================================================

/// <summary>
///     Custom action for uploading V2000 molfile (SDF).
/// </summary>
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

/// <summary>
///     Upload file button for uploading V2000 molfile (SDF).
/// </summary>
let private uploadFileButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                prop.children [
                    Html.span "Select SDF (Mol V2000) file"
                    uploadFileEvent dispatch
                ]
            ]
        ]
    ]

// ============================================================================
// GUI element: reset viewer button.
// ============================================================================

/// <summary>
///     Reset viewer button.
/// </summary>
let private resetViewerButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                prop.children [
                    Html.span "Reset viewer"
                ]
                prop.onClick (fun _ -> ResetViewer |> dispatch)
            ]
        ]
    ]

// ============================================================================
// GUI element: download button.
// ============================================================================

/// <summary>
///     Reset viewer button.
/// </summary>
let private downloadButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                prop.children [
                    Html.span "Download SVG"
                ]
                prop.onClick (fun _ -> DownloadSvg |> dispatch)
            ]
        ]
    ]

// ============================================================================
// GUI element: show hydrogens button.
// ============================================================================

/// <summary>
///     Reset viewer button.
/// </summary>
let private showHydrogensButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                prop.children [
                    Html.span "Toggle show hydrogens"
                ]
                prop.onClick (fun _ -> ToggleShowHydrogenAtoms |> dispatch)
            ]
        ]
    ]

// ============================================================================
// GUI element: change depiction button.
// ============================================================================

/// <summary>
///     Reset viewer button.
/// </summary>
let private changeDepictionButton dispatch =
    Html.div [
        prop.children [
            Bulma.button.a [
                prop.children [
                    Html.span "Toggle depiction"
                ]
                prop.onClick (fun _ -> ToggleDepiction |> dispatch)
            ]
        ]
    ]

// ============================================================================
// GUI element: SVG viewer.
// ============================================================================

/// <summary>
///     SVG viewer.
/// </summary>
let private svgViewer (dispatch: Msg -> unit) model =
    let svg =
        match model.EncodedSvg with
        | e when e.Length = 0 -> ""
        | _ -> $"data:image/svg+xml;base64,{model.EncodedSvg}"

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

// ============================================================================
// Main view.
// ============================================================================

/// <container>
///     Main view.
/// </container>
let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "cinemol"
        prop.children [
            Html.div [
                resetViewerButton dispatch
                uploadFileButton dispatch
                downloadButton dispatch
                showHydrogensButton dispatch
                changeDepictionButton dispatch
            ]
            Html.div [
                svgViewer dispatch model
            ]
        ]
    ]