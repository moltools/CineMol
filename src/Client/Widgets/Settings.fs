module Client.Widgets.Settings

open Feliz
open Feliz.Bulma
open Elmish
open Fable.React


type Model = float option

type Msg =
    | SetShowHydrogenAtoms
    | SetRotation of float
    | Oops of exn

type ExternalMessage =
    | NoOp
    | GotShowHydrogenAtoms
    | GotRotation of float

let init () =
    let model = None
    model, Cmd.none

let update (msg : Msg) (model : Model) : Model * ExternalMessage =
    match msg with
    | SetShowHydrogenAtoms ->
        model, ExternalMessage.GotShowHydrogenAtoms
    | SetRotation rotation ->
        model, (ExternalMessage.GotRotation rotation)
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
        Html.div [
            Slider.slider [
                slider.isFullWidth
                slider.isCircle
                slider.isLarge
                color.isBlack
                prop.onChange (fun (ev: Browser.Types.Event) -> (SetRotation (float ev.Value) |> dispatch))
            ]
        ]
    ]