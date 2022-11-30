module Client.Widgets.Settings

open Feliz
open Feliz.Bulma
open Elmish
open Fable.React


type Model = float option

type Msg =
    | SetShowHydrogenAtoms
    | SetDepiction
    | SetXRotation of float
    | SetYRotation of float
    | SetZRotation of float
    | Oops of exn

type ExternalMessage =
    | NoOp
    | GotShowHydrogenAtoms
    | GotDepiction
    | GotXRotation of float
    | GotYRotation of float
    | GotZRotation of float

let init () =
    let model = None
    model, Cmd.none

let update (msg : Msg) (model : Model) : Model * ExternalMessage =
    match msg with
    | SetShowHydrogenAtoms ->
        model, ExternalMessage.GotShowHydrogenAtoms
    | SetDepiction ->
        model, ExternalMessage.GotDepiction
    | SetXRotation rotation ->
        model, (ExternalMessage.GotXRotation rotation)
    | SetYRotation rotation ->
        model, (ExternalMessage.GotYRotation rotation)
    | SetZRotation rotation ->
        model, (ExternalMessage.GotZRotation rotation)
    | Msg.Oops _ -> model, ExternalMessage.NoOp  // ignored for now

let view (model : Model) (dispatch : Msg -> unit) =
    Bulma.content [
        Html.div [
            Switch.checkbox [
                prop.id "hydrogen-switch"
                prop.onChange (fun (_ : Browser.Types.Event) -> SetShowHydrogenAtoms |> dispatch)
                color.isSuccess
            ]
            Html.label [
                prop.htmlFor "hydrogen-switch"
                prop.text "Show hydrogen atoms"
            ]
        ]
//        Html.div [
//            Switch.checkbox [
//                prop.id "depiction-switch"
//                prop.onChange (fun (_ : Browser.Types.Event) -> SetDepiction |> dispatch)
//                color.isSuccess
//            ]
//            Html.label [
//                prop.htmlFor "depiction-switch"
//                prop.text "Change depiction"
//            ]
//        ]
        Html.div [
            Slider.slider [
                slider.isFullWidth
                slider.isCircle
                slider.isLarge
                color.isBlack
                prop.onChange (fun (ev: Browser.Types.Event) -> (SetXRotation (float ev.Value) |> dispatch))
            ]
        ]
        Html.div [
            Slider.slider [
                slider.isFullWidth
                slider.isCircle
                slider.isLarge
                color.isBlack
                prop.onChange (fun (ev: Browser.Types.Event) -> (SetYRotation (float ev.Value) |> dispatch))
            ]
        ]
        Html.div [
            Slider.slider [
                slider.isFullWidth
                slider.isCircle
                slider.isLarge
                color.isBlack
                prop.onChange (fun (ev: Browser.Types.Event) -> (SetZRotation (float ev.Value) |> dispatch))
            ]
        ]
    ]