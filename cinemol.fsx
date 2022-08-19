// Author: David Meijer
// Description: 3D SVG drawer for small molecules.
// Usage: dotnet fsi cinemol.fsx <path to SDF/Mol V2000 file>
open System
open System.IO
open System.Text.RegularExpressions

exception InputError of string 

let args : string array = fsi.CommandLineArgs |> Array.tail
let input_path : string = args.[0]

let atomLine : string = @"^\s{1,}([-+]?[0-9]*\.?[0-9]+)\s{1,}([-+]?[0-9]*\.?[0-9]+)\s{1,}([-+]?[0-9]*\.?[0-9]+)\s{1,}(\w+)\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+\s{1,}[-+]?[0-9]*\.?[0-9]+$"

let (|AtomLine|_|) input =
    let m = Regex.Match(input, atomLine)
    if m.Success then Some(List.tail [ for g in m.Groups -> g.Value ])
    else None 

type Atom = 
    | C | N | O | S | H
    member this.Radius : float =
        match this with 
        | H -> 1.0
        | _ -> 1.2

    member this.Color : string =
        match this with 
        | C -> "#808080"
        | N -> "#0000FF"
        | O -> "#FF0000"
        | S -> "#FFFF00"
        | H -> "#D3D3D3"

type Coords = { 
    X : float 
    Y : float
    Z : float 
}

let identifyAtom (atom : string) : Atom =
    match atom with 
    | "C" -> C 
    | "N" -> N
    | "O" -> O
    | "S" -> S
    | "H" -> H
    | _ -> raise <| InputError($"unknown atom {atom}")

let tryParseFloat (s : string) : float option =
    try s |> float |> Some 
    with :? FormatException -> None 

let castToFloat (s : string) : float =
    match tryParseFloat s with 
    | Some f -> f 
    | None -> raise <| InputError("coordinate is not a float")

let parse_sdf (sdf_path : string) : (int * Atom * Coords) array =
    use reader = new IO.StreamReader(sdf_path) 
    let lines = [| while not reader.EndOfStream do yield reader.ReadLine() |]
    
    let moleculeCount =
        lines 
        |> Array.map (fun l -> l.Contains("$$$$") = true) 
        |> Array.filter id
        |> Array.length

    match moleculeCount with
    | x when x <> 1 -> raise <| InputError("multiple molecules in input file")
    | _ -> ()

    let mutable atomCount = 0
    [|
        for l in lines do
            match l with 
            | AtomLine [ x; y; z; symbol] -> 
                atomCount <- atomCount + 1
                (
                    atomCount,
                    identifyAtom symbol, 
                    {
                        X = castToFloat x
                        Y = castToFloat y
                        Z = castToFloat z
                    }
                )
            | _ -> ()
    |]

let writeAtomStyle ((idx, atom, _) : int * Atom * Coords) : string =
    $".atom-{idx}{{fill:{atom.Color};}}"

let writeAtom ((idx, atom, coords) : int * Atom * Coords) : string =
    $"<circle class=\"atom-{idx}\" cx=\"{coords.X}\" cy=\"{coords.Y}\" r=\"{atom.Radius}\"/>"

let main =
    let camera = { X = 0.0; Y = 0.0; Z = 10.0 }

    let atoms = 
        parse_sdf input_path
        |> Array.filter (fun ((idx, atom, coords) : int * Atom * Coords ) -> atom <> H)

    let xMin = -10
    let yMin = -10
    let width = 20
    let height = 20
    
    use w = new StreamWriter("out.svg")

    w.WriteLine("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    w.WriteLine($"<svg id=\"Layer_1\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{xMin} {yMin} {width} {height}\">")
    w.WriteLine("<defs>")
    w.WriteLine("<style>")

    atoms 
    |> Array.map writeAtomStyle 
    |> Array.map (fun l -> w.WriteLine(l))
    |> ignore

    w.WriteLine("</style>")
    w.WriteLine("</defs>")

    atoms
    |> Array.map writeAtom
    |> Array.map (fun l -> w.WriteLine(l))
    |> ignore 
    
    w.WriteLine("</svg>")

    0

[<EntryPoint>]
main |> printfn "%A"
