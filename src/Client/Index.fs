module Index

open Browser
open Browser.Types
open Elmish
open Fable
open Fable.React
open Fable.React.Props
open Fable.Remoting.Client
open Fable.Core
open Feliz
open Feliz.Bulma
open Fulma
open Shared

let encodeURI (svg : string) : string =
//    JS.encodeURIComponent svg
    svg

type Model =
    { Todos: Todo list
      Input: string
      Sdf: string
      Svg: string
      Encoded: string }

type Msg =
    | GotTodos of Todo list
    | SetInput of string
    | AddTodo
    | AddedTodo of Todo
    | UploadSdf of name : string * content : string
    | Render
    | GotEncoding of string

let todosApi =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.buildProxy<ITodosApi>

let init () : Model * Cmd<Msg> =
    let model =
        { Todos = []
          Input = ""
          Sdf = ""
          Svg = ""
          Encoded = "" }
//    let cmd = Cmd.OfAsync.perform todosApi.getTodos () GotTodos
//    let cmd = Cmd.OfAsync.perform todosApi.render (encodeURI model.Svg) GotEncoding
    let cmd = Cmd.none
    model, cmd

let update (msg: Msg) (model: Model) : Model * Cmd<Msg> =
    match msg with
    | GotTodos todos -> { model with Todos = todos }, Cmd.none
    | SetInput value -> { model with Input = value }, Cmd.none
    | AddTodo ->
        let todo = Todo.create model.Input
        let cmd = Cmd.OfAsync.perform todosApi.addTodo todo AddedTodo
        { model with Input = "" }, cmd
    | AddedTodo todo -> { model with Todos = model.Todos @ [ todo ] }, Cmd.none
    | UploadSdf (_, content) -> { model with Sdf = content }, Cmd.ofMsg Render
    | Render -> model, Cmd.OfAsync.perform todosApi.render (encodeURI model.Sdf) GotEncoding
    | GotEncoding encoded -> { model with Encoded = encoded }, Cmd.none

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
            )
        ]
    ]

let private uploadFileButton dispatch =
    Html.div [
        prop.className "action-button"
        prop.children [
            Bulma.button.a [
                button.isOutlined
                button.isFullWidth
                prop.children [
                    Html.span "select SDF/Mol V2000 file"
                    uploadFileEvent dispatch
                ]
            ]
        ]
    ]

let private svgViewer model dispatch =
    match model.Encoded with
    | e when e.Length = 0 -> Html.div [ prop.className "viewer" ]
    | _ ->
        Html.div [
            prop.className "viewer"
            prop.children [
                img [ Src $"data:image/svg+xml;base64,{model.Encoded}"]
            ]
        ]

let view (model: Model) (dispatch: Msg -> unit) =
    Html.div [
        prop.className "cinemol"
        prop.children [
            uploadFileButton dispatch
            svgViewer model dispatch
        ]
    ]