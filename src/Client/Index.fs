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
type ViewerBackgroundStyle = | Dark | Light
with member x.toHex = match x with | Dark -> "#000000" | Light -> "#FFFFFF"

type Model = {
    Sdf: string
    Svg: string
    EncodedSvg: string
    ViewBox: ViewBox option
    DrawOptions: DrawOptions
    Rotation: Rotation
    Zoom: Zoom
    DragTarget: DragTarget
    ViewerBackgroundStyle: ViewerBackgroundStyle
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
            ViewerBackgroundStyle = Dark
        }

type Msg =
    /// GUI buttons.
    | UploadSdf of name: string * content: string
    | ResetViewer
    | DownloadSvg
    | ToggleShowHydrogenAtoms
    | ToggleDepiction
    | ToggleBackgroundStyle

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
        let newModel = { Model.init with ViewerBackgroundStyle = model.ViewerBackgroundStyle }
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
            elif model.DrawOptions.Depiction = Depiction.BallAndStick then
                Depiction.Wire
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

    /// Toggle backgroundc color
    | ToggleBackgroundStyle ->
        let toggle = if model.ViewerBackgroundStyle = Dark then Light else Dark
        let newModel = { model with ViewerBackgroundStyle = toggle }
        newModel, Cmd.none

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

    | MouseDragStarted _ ->
        { model with DragTarget = Dragging }, Cmd.none

    | MouseDragEnded ->
        { model with DragTarget = NoTarget }, Cmd.none

    | MouseDrag (position: MousePosition) ->
        match model.DragTarget with
        | Dragging -> model, Cmd.ofMsg (SetRotation position)
        | _ -> model, Cmd.none

    /// Mouse scroll message for rendered, encoded SVG in model.
    | WheelScroll (wheel: WheelPosition) ->
        let incr = if wheel.Delta < 0.0 then -0.05 else +0.05
        let newRatio =
            let newRatio = model.Zoom.Ratio + incr
            if newRatio < 0.0 then 0.0 else newRatio
        let newModel = { model with Zoom = { model.Zoom with Ratio = newRatio } }
        newModel, Cmd.OfAsync.perform render newModel.renderArgs GotEncoding

// ============================================================================
// Buttons
// ============================================================================
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

                    /// Setting value to empty string removes cached filename
                    /// again and makes sure you can upload the same file
                    /// (after refresh) for a second time.
                    (ev.target :?> HTMLInputElement).value <- ""
            )
        ]
    ]

let private uploadFileButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.children [
            Html.i [ prop.className "fas fa-upload" ]
            uploadFileEvent dispatch
        ]
    ]

let private resetViewerButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.children [ Html.i [ prop.className "fas fa-sync" ] ]
        prop.onClick (fun _ -> ResetViewer |> dispatch)
    ]

let private downloadButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.children [ Html.i [ prop.className "fas fa-download" ] ]
        prop.onClick (fun _ -> DownloadSvg |> dispatch)
    ]

let private showHydrogensButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.style [ style.fontWeight 700; style.color "#000000" ]
        prop.children [ Html.span "Hs" ]
        prop.onClick (fun _ -> ToggleShowHydrogenAtoms |> dispatch)
    ]

let private changeDepictionButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.children [ Html.i [ prop.className "fas fa-eye" ] ]
        prop.onClick (fun _ -> ToggleDepiction |> dispatch)
    ]

let private changeBackgroundStyleButton dispatch =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.children [ Html.i [ prop.className "fas fa-adjust" ] ]
        prop.onClick (fun _ -> ToggleBackgroundStyle |> dispatch)
    ]

let private reportBugButton =
    Bulma.button.a [
        prop.className "sidebar-button"
        prop.href "https://github.com/moltools/cinemol/issues"
        prop.rel "noreffer noopener"
        prop.target "_blank"
        prop.children [
            Html.i [ prop.className "fas fa-bug" ]
        ]
    ]


// ============================================================================
// GUI element: SVG viewer.
// ============================================================================
let private svgViewer (dispatch: Msg -> unit) model =
    let svg =
        match model.EncodedSvg with
        | e when e.Length = 0 ->
            let emptySvg = $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg id=\"Layer_1\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{0} {0} {100} {100}\"><style></style></svg>"
            let encoded = emptySvg |> toBase64String
            $"data:image/svg+xml;base64,{encoded}"
        | _ ->
            $"data:image/svg+xml;base64,{model.EncodedSvg}"

    let size =
        let width = int window.innerWidth - 50
        let height = int window.innerHeight
        if width < height then
            [ style.width width; style.height width ]
        else
            [ style.width height; style.height height ]

    Html.div [
        prop.className "viewer-window"
        prop.style size
        prop.onMouseDown (fun ev ->
            ev.preventDefault()
            let coordsMouseDown = { X = ev.pageX; Y = ev.pageY }
            dispatch (MouseDragStarted coordsMouseDown))
        prop.children [
            img [
                Class "viewer-svg"
                Style [ BackgroundColor model.ViewerBackgroundStyle.toHex ]
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
        prop.style [ style.backgroundColor model.ViewerBackgroundStyle.toHex ]
        prop.children [
            Html.div [
                prop.className "sidebar"
                prop.children [
                    resetViewerButton dispatch
                    uploadFileButton dispatch
                    downloadButton dispatch
                    showHydrogensButton dispatch
                    changeDepictionButton dispatch
                    changeBackgroundStyleButton dispatch
                    reportBugButton
                ]
            ]
            svgViewer dispatch model
        ]
    ]