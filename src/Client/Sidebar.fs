module Client.Sidebar

open Client.Widgets.General
open Elmish
open Feliz
open Feliz.Bulma
open Fable.React


type Html with
    static member inline ofOption (element : ReactElement option) : ReactElement =
        match element with
        | Some element -> element
        | None -> Html.none

type Model =
    { IsExpanded : bool
      WidgetState : Set<string>
      Settings : Widgets.Settings.Model
      General : Widgets.General.Model }

type Msg =
    | GeneralMsg of Widgets.General.Msg
    | SettingsMsg of Widgets.Settings.Msg
    | ToggleWidget of string
    | ToggleState

type ExternalMsg =
    | Reset
    | UploadSdf of name : string * src : string
    | DownloadSvg
    | GotShowHydrogenAtoms
    | GotDepiction
    | NoOp

let init () =
    let settingsModel, settingsCmd = Widgets.Settings.init()
    let cmd = Cmd.map SettingsMsg settingsCmd
    let model =
        { IsExpanded = true
          WidgetState = Set.empty.Add("General").Add("Settings")
          General = Widgets.General.init()
          Settings = settingsModel }
    model, cmd

let update (msg : Msg) (model : Model) : Model * Cmd<Msg> * ExternalMsg =
    match msg with
    | GeneralMsg msg ->
        let generalModel, externalMsg = update msg model.General
        let externalMsg =
            match externalMsg with
            | ExternalMessage.NoOp -> NoOp
            | ExternalMessage.Reset -> Reset
            | ExternalMessage.UploadSdf (name, src) -> UploadSdf (name, src)
            | ExternalMessage.DownloadSvg -> DownloadSvg
        { model with General = generalModel }, Cmd.none, externalMsg
    | SettingsMsg msg ->
        let settingsModel, externalMsg = Widgets.Settings.update msg model.Settings
        let externalMsg =
            match externalMsg with
            | Widgets.Settings.NoOp -> ExternalMsg.NoOp
            | Widgets.Settings.GotShowHydrogenAtoms -> ExternalMsg.GotShowHydrogenAtoms
            | Widgets.Settings.GotDepiction -> ExternalMsg.GotDepiction
        { model with Settings = settingsModel }, Cmd.map SettingsMsg Cmd.none, externalMsg
    | ToggleWidget id ->
        let newWidgetState =
            if model.WidgetState.Contains id then model.WidgetState.Remove id
            else model.WidgetState.Add id
        { model with WidgetState = newWidgetState }, Cmd.none, NoOp
    | ToggleState -> { model with IsExpanded = not model.IsExpanded }, Cmd.none, NoOp

let private renderExpandedWidgets (states : Set<string>) dispatch (title, icon, widget : ReactElement, maxHeight : int option) =
    let baseView headerIcon content =
        Bulma.card [
            Bulma.cardHeader [
                prop.onClick (fun _ -> ToggleWidget title |> dispatch )
                prop.children [
                    Bulma.cardHeaderTitle.div [
                        Bulma.icon [
                            prop.style [
                                Feliz.style.marginRight (length.em 0.5)
                            ]
                            prop.children [
                                icon
                            ]
                        ]
                        Html.text title
                    ]
                    Bulma.cardHeaderIcon.span [
                        Bulma.icon [
                            icon
                        ]
                    ]
                ]
            ]
            Html.ofOption content
        ]
    if states.Contains title then
            baseView (Html.i [ prop.className "fas fa-angle-down" ]) None
        else
            let props : IReactProperty list =
                match maxHeight with
                | Some maxHeight ->
                    [
                        prop.style [
                            Feliz.style.maxHeight (length.px maxHeight)
                            Feliz.style.overflowY.auto
                        ]
                    ]
                | None -> []
            baseView (Html.i [ prop.className "fas fa-angle-up" ]) (Some (Bulma.cardContent [yield! props; prop.children widget]))

let renderCollapsedWidgets dispatch (title, faIcon, widget: ReactElement, maxHeight: int option) =
    let props =
        match maxHeight with
        | Some maxHeight ->
            [
                prop.style [
                    Feliz.style.maxHeight maxHeight
                    Feliz.style.overflowY.auto
                ]
            ]
        | None -> []

    Html.div [
        prop.className "item"
        prop.children [
            Bulma.icon [
                prop.className "is-large"
                prop.children [
                    faIcon
                ]
            ]
            Bulma.card [
                prop.className "item-content"
                prop.children [
                    Bulma.cardHeader [
                        prop.onClick (fun _ -> ToggleWidget title |> dispatch )
                        prop.children [
                            Bulma.cardHeaderTitle.div title
                        ]
                    ]
                    Bulma.cardContent [
                        yield! props
                        prop.children widget
                    ]
                ]
            ]
        ]
    ]


let private renderWidgets model dispatch (title, icon, widget, maxHeight) =
    match model.IsExpanded with
    | true ->
        renderExpandedWidgets model.WidgetState dispatch (title, icon, widget, maxHeight)
    | false ->
        renderCollapsedWidgets dispatch (title, icon, widget, maxHeight)

let private collapseButton dispatch =
    Bulma.card [
        prop.onClick (fun _ -> dispatch ToggleState)
        prop.children [
            Bulma.cardHeader [
                Bulma.cardHeaderTitle.div "Collapse"
                Bulma.cardHeaderIcon.span [
                    Bulma.icon [
                        (Html.i [ prop.className "fas fa-angle-double-left" ])
                    ]
                ]
            ]
        ]
    ]


let private expandButton dispatch =
    Bulma.card [
        prop.onClick (fun _ -> dispatch ToggleState)
        prop.children [
            Bulma.cardHeader [
                Bulma.cardHeaderIcon.span [
                    Bulma.icon [
                        (Html.i [ prop.className "fas fa-angle-double-right" ])
                    ]
                ]
            ]
        ]
    ]

let private sidebarContainer dispatch (sections : ReactElement list) =
    Html.div [
        prop.className "sidebar is-expanded"
        prop.children [
            Html.div [
                prop.className "brand"
                prop.children [
                    Bulma.title.h4 ""
                ]
            ]
            Html.div [
                prop.className "widgets-list"
                prop.children sections
            ]
//            collapseButton dispatch
        ]
    ]

let view (model : Model) (dispatch : Msg -> unit) : ReactElement =
    let widgets =
        [
            if model.IsExpanded then
                "General", (Html.i [ prop.className "fas fa-book" ]), Widgets.General.viewExpanded model.General (GeneralMsg >> dispatch), None
            "Settings", (Html.i [ prop.className "fas fa-cog" ]), Widgets.Settings.view model.Settings (SettingsMsg >> dispatch), Some 300
        ]
        |> List.map (renderWidgets model dispatch)

    if model.IsExpanded then
        sidebarContainer dispatch widgets
    else
        let generalCollapsedView =
            viewCollapsed model.General (GeneralMsg >> dispatch)

        Html.div [
            prop.className "sidebar is-collapsed"
            prop.children [
                viewModalResetConfirmation model.General (GeneralMsg >> dispatch)

                Html.div [
                    prop.className "brand"
                ]
                Html.div [
                    prop.className "widgets-list"
                    prop.children (generalCollapsedView::widgets)
                ]
                expandButton dispatch
            ]
        ]
