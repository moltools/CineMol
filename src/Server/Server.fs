module Server

open Fable.Remoting.Server
open Fable.Remoting.Giraffe
open Saturn

open Shared
open Cinemole.Helpers
open Cinemole.Types
open Cinemole.Drawing

let cinemoleApi =
    { render = fun assignment -> async {
        let depiction =
            match assignment.Settings.Depiction with
            | Shared.Filled -> Filled
            | Shared.BallAndStick -> BallAndStick

        let svg, viewBox : string * ViewBox =
            assignment.Sdf
            |> draw assignment.Settings.ViewBox
                    depiction
                    assignment.Settings.ShowHydrogenAtoms
                    assignment.Settings.XRotation
                    assignment.Settings.YRotation
                    assignment.Settings.ZRotation

        let encodedSvg : string = svg |> toBase64String

        return (svg, encodedSvg, viewBox) } }

let webApp =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.fromValue cinemoleApi
    |> Remoting.buildHttpHandler

let app =
    application {
        url "http://*:8085"
        use_router webApp
        memory_cache
        use_static "public"
        use_gzip }

[<EntryPoint>]
let main _ =
    run app
    0