module Client.CineMol.Parsing

open System
open System.Text.RegularExpressions

open Types

/// <summary>
///     Exception to throw when parser encounters error when parsing.
/// </summary>
exception ParserError of string

/// <summary>
///     Constructs regular expression for atom line in V2000 molfile.
/// </summary>
/// <returns>
///     Regular expression for atom line in V2000 molfile that matches four
///     groups (in the following order):
///         1) X coordinate
///         2) Y coordinate
///         3) Z coordinate
///         4) atom type
/// </returns>
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

/// <summary>
///     Active pattern for matching and parsing atom line in V2000 molfile
/// </summary>
/// <param name="line">
///     Input line to apply atom line regular expression on.
/// </param>
/// <returns>
///     List of matched groups, if there was a match with the regular
///     expression.
/// </returns>
let (|AtomLine|_|) (line: string) : string list option =
    let m = Regex.Match(line, atomLine)
    if m.Success then
        Some(List.tail [ for g in m.Groups -> g.Value ])
    else
        None

/// <summary>
///     Cast atom symbol to Atom.
/// </summary>
/// <param name="atom">
///     Atom symbol.
/// </param>
/// <returns>
///     Casted atom type.
/// </returns>
let tryCastToAtom (atom: string) : Atom =
    match atom with
    | "C" -> C | "N" -> N | "O" -> O | "S" -> S | "H" -> H
    | _ -> Unknown

/// <summary>
///     Try to cast float string to float.
/// </summary>
/// <param name="s">
///     Float string.
/// </param>
/// <returns>
///     Casted float, if successfully cast. Otherwise, float 0.0.
/// </returns>
let tryCastToFloat (s: string) : float =
    try s |> float
    with :? FormatException -> 0.0

/// <summary>
///     Parse molecules from V2000 molfile.
/// </summary>
/// <param name="sdf">
///     SDF string to parse molecules from.
/// </param>
/// <returns>
///     Parsed molecules from SDF string.
/// </returns>
let parseSdf (sdf: string) : Molecule[] =
    let mutable atoms: AtomInfo list = []
    let mutable atomCount: int = 0
    [|
        for line in sdf.Split [|'\n'|] do
            match line with
            | line when line.Contains("$$$$") = true ->
                yield { Atoms = atoms |> List.toArray }
                atomCount <- 0
                atoms <- []
            | AtomLine [ x; y; z; symbol] ->
                atomCount <- atomCount + 1
                let atomType = tryCastToAtom symbol
                let center: Point = {
                    X = tryCastToFloat x
                    Y = tryCastToFloat y
                    Z = tryCastToFloat z
                }
                let radius: Radius = atomType.Radius
                let atom = createAtom atomCount atomType center radius
                atoms <- atoms @ [ atom ]
            | _ -> ()
    |]
