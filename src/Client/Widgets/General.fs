module Client.Widgets.General

open Browser.Types
open Fable.React.Props
open Feliz
open Feliz.Bulma
open Browser
open Fulma

type ResetState =
    | Default
    | Confirm

type Model =
    { ResetState : ResetState }

type Msg =
    | AskReset
    | ConfirmReset
    | CancelReset
    | UploadSdf of name : string * src : string
    | DownloadScript

[<RequireQualifiedAccess>]
type ExternalMessage =
    | NoOp
    | Reset
    | UploadSdf of name : string * src : string
    | DownloadSvg

type ActionButtonType =
    | Default of Msg
    | RequestFileInput

let init () = { ResetState = ResetState.Default }

let update (msg : Msg) (model : Model) : Model * ExternalMessage =
    match msg with
    | AskReset -> { model with ResetState = Confirm }, ExternalMessage.NoOp
    | ConfirmReset -> { model with ResetState = ResetState.Default }, ExternalMessage.Reset
    | CancelReset -> { model with ResetState = ResetState.Default }, ExternalMessage.NoOp
    | UploadSdf (name, src) -> model, ExternalMessage.UploadSdf (name, src)
    | DownloadScript -> model, ExternalMessage.DownloadSvg

let private fileInputEvent dispatch =
    Fulma.File.input [
        Props [
            Style [
                Cursor "pointer"
                Display DisplayOptions.Inherit
            ]
            OnInput (fun ev ->
                let file = (ev.target :?> HTMLInputElement).files.[0]
                let name = file.name
                let reader = FileReader.Create()
                reader.onload <- fun _ ->
                    let content = reader.result :?> string
                    (name, content) |> UploadSdf |> dispatch
                reader.readAsText file)
        ]
    ]

let viewCollapsed (model : Model) dispatch =
    let uploadIcon = [ Html.i [ prop.className "fas fa-star" ] ]
    let downloadIcon = [ Html.i [ prop.className "fas fa-star" ] ]
    let refreshIcon = [ Html.i [ prop.className "fas fa-star" ] ]

    let actionButton (wrapped : ActionButtonType) (faIcon : ReactElement list) =
        match wrapped with
        | Default msg ->
            Html.div [
                prop.className "action-button"
                prop.children [
                    Bulma.button.a [
                        button.isOutlined
                        prop.onClick (fun _ -> dispatch msg)
                        prop.children [
                            Bulma.icon [
                                icon.isLarge
                                prop.children faIcon
                            ]
                        ]
                    ]
                ]
            ]
        | RequestFileInput ->
            Html.div [
                prop.className "action-button"
                prop.children [
                    Bulma.button.a [
                        prop.style [
                            Feliz.style.width 35
                        ]
                        button.isOutlined
                        prop.children [
                            Bulma.icon [
                                prop.style [
                                    Feliz.style.marginRight -8
                                ]
                                icon.isLarge
                                prop.children faIcon
                            ]
                            fileInputEvent dispatch
                        ]
                    ]
                ]
            ]

    Html.div [
        prop.className "actions-area"
        prop.children [
            actionButton (Default AskReset) refreshIcon
            actionButton RequestFileInput uploadIcon
            actionButton (Default DownloadScript) downloadIcon
        ]
    ]

let viewExpanded (model : Model) dispatch =
    let renderItem (wrapped : ActionButtonType) (text : string) =

        match wrapped with
        | Default msg ->
            Bulma.field.div [
                field.hasAddons
                prop.children [
//                    Bulma.control.div [
//                        Bulma.button.a [
//                            prop.onClick (fun _ -> dispatch msg)
//                        ]
//                    ]
                    Bulma.control.div [
                        control.isExpanded
                        prop.children [
                            Bulma.button.a [
                                prop.onClick (fun _ -> dispatch msg)
                                button.isText
                                button.isFullWidth
                                prop.children [
                                    Html.span text
                                ]
                            ]
                        ]
                    ]
                ]
            ]

        | RequestFileInput ->
            Bulma.field.div [
                field.hasAddons
                prop.children [
//                    Bulma.control.div [
//                        Bulma.button.a [
//                            prop.children [
//                                fileInputEvent dispatch
//                            ]
//                        ]
//                    ]
                    Bulma.control.div [
                        control.isExpanded
                        prop.children [
                            Bulma.button.a [
                                button.isText
                                button.isFullWidth
                                prop.children [
                                    Html.span text
                                    fileInputEvent dispatch
                                ]
                            ]
                        ]
                    ]
                ]
            ]

    let content =
        match model.ResetState with
        | ResetState.Default ->
            Html.div [
                renderItem (Default AskReset) " Refresh renderer"
                renderItem RequestFileInput "Upload SDF (V2000)"
                renderItem (Default DownloadScript) "Download SVG"
            ]


        | Confirm ->
            Bulma.field.div [
                Bulma.help [
                    color.isBlack
                    prop.text "Please, confirm to reset"
                ]

                Bulma.field.div [
                    field.hasAddons
                    prop.children [
                        Bulma.control.div [
                            Bulma.button.a [
                                prop.onClick (fun _ -> dispatch ConfirmReset)
                                color.isSuccess
                                prop.children [
                                    Bulma.icon [
                                        Html.i [ prop.className "fas fa-check" ]
                                    ]
                                    Html.span "Confirm"
                                ]
                            ]
                        ]

                        Bulma.control.div [
                            Bulma.button.a [
                                prop.onClick (fun _ -> dispatch CancelReset)
                                color.isDanger
                                prop.children [
                                    Bulma.icon [
                                        Html.i [ prop.className "fas fa-times" ]
                                    ]
                                    Html.span "Cancel"
                                ]
                            ]
                        ]
                    ]
                ]
            ]

    Bulma.content content

let viewModalResetConfirmation (model: Model) dispatch =
    Bulma.modal [
        if (model.ResetState = Confirm) then
            modal.isActive
        prop.style [
            // Make sure to be on top of everything.
            style.zIndex 20000
        ]
        prop.children [
            Bulma.modalBackground [
                prop.onClick (fun _ -> dispatch CancelReset)
            ]

            Bulma.modalContent [
                Html.div [
                    prop.className "reset-confirmation-modal"
                    prop.children [
                        Html.div [
                            prop.className "reset-confirmation-modal-content"
                            prop.children [
                                Html.span [
                                    prop.className "reset-confirmation-modal-content-text"
                                    prop.text "Please, confirm to reset"
                                ]

                                Html.div [
                                    prop.className "reset-confirmation-modal-content-foot"
                                    prop.children [
                                        Bulma.field.div [
                                            prop.children [
                                                Bulma.field.div [
                                                    field.hasAddons
                                                    prop.children [
                                                        Bulma.control.div [
                                                            Bulma.button.a [
                                                                prop.onClick (fun _ -> dispatch ConfirmReset)
                                                                color.isSuccess
                                                                prop.children [
                                                                    Bulma.icon [
                                                                        (Html.i [ prop.className "fas fa-star" ])
                                                                    ]

                                                                    Html.span "Confirm"
                                                                ]
                                                            ]
                                                        ]

                                                        Bulma.control.div [
                                                            Bulma.button.a [
                                                                prop.onClick (fun _ -> dispatch CancelReset)
                                                                color.isDanger
                                                                prop.children [
                                                                    Bulma.icon [
                                                                        (Html.i [ prop.className "fas fa-star" ])
                                                                    ]

                                                                    Html.span "Cancel"
                                                                ]
                                                            ]
                                                        ]
                                                    ]
                                                ]
                                            ]
                                        ]
                                    ]
                                ]
                            ]
                        ]
                    ]
                ]
            ]
        ]
    ]