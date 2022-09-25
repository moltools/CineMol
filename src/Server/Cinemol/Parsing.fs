namespace Cinemol

module Parsing =

    open System
    open System.Text.RegularExpressions

    open Cinemol.Types
    open Cinemol.ErrorHandling

    let atomLine : string =
        let s = @"\s{1,}"
        let d = @"[-+]?[0-9]*\.?[0-9]+"
        let d_cap = $"({d})"
        let w_cap = $@"(\w+)"
        [ "$" ]
        |> (@) [ for _ in [ 0 .. 11 ] do yield s + d ]
        |> (@) [ s + w_cap]
        |> (@) [ for _ in [ 0 .. 2 ] do yield s + d_cap ]
        |> (@) [ "^" ]
        |> String.concat ""

    let (|AtomLine|_|) input =
        let m = Regex.Match(input, atomLine)
        if m.Success then Some(List.tail [ for g in m.Groups -> g.Value ]) else None

    let identifyAtom (atom: string) : Atom =
        match atom with
        | "C" -> C | "N" -> N | "O" -> O | "S" -> S | "H" -> H
        | _ -> raise <| InputError($"unknown atom {atom}")

    let tryParseFloat (s: string) : float option =
        try s |> float |> Some
        with :? FormatException -> None

    let castToFloat (s: string) : float =
        match tryParseFloat s with | Some f -> f | None -> raise <| InputError("coordinate is not a float")

    let parseSdf (sdf: string) : AtomInfo array =
        let lines = sdf.Split [|'\n'|]
        let moleculeCount = lines |> Array.map (fun l -> l.Contains("$$$$") = true) |> Array.filter id |> Array.length
        match moleculeCount with | x when x > 1 -> raise <| InputError("multiple molecules in input file") | _ -> ()

        let mutable atomCount = 0
        [| for l in lines do
            match l with
            | AtomLine [ x; y; z; symbol] ->
                atomCount <- atomCount + 1
                let atomType = identifyAtom symbol
                let center: Point = { X = castToFloat x; Y = castToFloat y; Z = castToFloat z }
                let radius: Radius = atomType.Radius
                createAtom atomCount atomType center radius
            | _ -> () |]