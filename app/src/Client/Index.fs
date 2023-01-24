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

type Model = { ID: string option }
    with
    static member init = { ID = None }

type Msg = | DoNothing

let init () : Model * Cmd<Msg> =
    Model.init, Cmd.none

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | DoNothing -> model, Cmd.none

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div []