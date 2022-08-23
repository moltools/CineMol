module Server

open Fable.Remoting.Server
open Fable.Remoting.Giraffe
open Saturn
open System
open System.Text.RegularExpressions

open Shared

// ============================================================================
// Error handling.
// ============================================================================
exception InputError of string

// ============================================================================
// Globals.
// ============================================================================
let dr1, dr2, dr3, dr4, dr5 = 0.0, 0.11, 0.34, 0.66, 1.0  // Diffusion rates

// ============================================================================
// Types.
// ============================================================================
type Gradient = Color * Color * Color * Color * Color

and Color =
    { R : int
      G : int
      B : int }
    member this.Diffuse (factor : float) : Color =
        { R = int ((float this.R) * factor)
          G = int ((float this.G) * factor)
          B = int ((float this.B) * factor) }
    member this.Gradient : Gradient =
        ( this.Diffuse (1.0 - dr1),
          this.Diffuse (1.0 - dr2),
          this.Diffuse (1.0 - dr3),
          this.Diffuse (1.0 - dr4),
          this.Diffuse (1.0 - dr5) )

type Index = int

type Atom =
    | C | N | O | S | H
    member this.Radius : float =
        match this with
        | H -> 1.0
        | _ -> 1.2
    member this.Color : Color =
        match this with
        | C -> { R = 128; G = 128; B = 128 }
        | N -> { R = 0; G = 0; B = 255 }
        | O -> { R = 255; G = 0; B = 0 }
        | S -> { R = 255; G = 255; B = 0 }
        | H -> { R = 220; G = 220; B = 220 }

type Coords =
    { X : float
      Y : float
      Z : float }
    static member (-) (c1 : Coords, c2 : Coords) : Coords =
        { X = c1.X - c2.X; Y = c1.Y - c2.Y; Z = c1.Z - c2.Z }
    static member Pow (c : Coords) (d : float) : Coords =
        { X = c.X ** d; Y = c.Y ** d; Z = c.Z ** d }
    static member Sum (c : Coords) : float =
        c.X + c.Y + c.Z

type AtomInfo = Index * Atom * Coords

type ViewBox = float * float * float * float

// ============================================================================
// Parsing input SDF/Mol V2000 file.
// ============================================================================
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
    if m.Success then Some(List.tail [ for g in m.Groups -> g.Value ])
    else None

let identifyAtom (atom : string) : Atom =
    match atom with
    | "C" -> C | "N" -> N | "O" -> O | "S" -> S | "H" -> H
    | _ -> raise <| InputError($"unknown atom {atom}")

let tryParseFloat (s : string) : float option =
    try s |> float |> Some
    with :? FormatException -> None

let castToFloat (s : string) : float =
    match tryParseFloat s with
    | Some f -> f
    | None -> raise <| InputError("coordinate is not a float")

let parse_sdf (sdf : string) : AtomInfo array =
    let lines = sdf.Split [|'\n'|]

    let moleculeCount =
        lines
        |> Array.map (fun l -> l.Contains("$$$$") = true)
        |> Array.filter id
        |> Array.length

    match moleculeCount with
    | x when x > 1 -> raise <| InputError("multiple molecules in input file")
    | _ -> ()

    let mutable atomCount = 0
    [| for l in lines do
        match l with
        | AtomLine [ x; y; z; symbol] ->
            atomCount <- atomCount + 1
            ( atomCount,
              identifyAtom symbol,
              { X = castToFloat x
                Y = castToFloat y
                Z = castToFloat z })
        | _ -> () |]

// ============================================================================
// Functions for manipulating objects.
// ============================================================================
let abs (v : float) : float =
    (v ** 2.0) ** 0.5

let dist (a : Coords) (b : Coords) : float =
    (Coords.Sum(Coords.Pow (a - b) 2.0)) ** 0.5

let rotateAxisX (c : Coords) (rad : float option) : Coords =
    match rad with
    | Some r ->
        { X = c.X
          Y = c.Y * Math.Cos(r) - c.Z * Math.Sin(r)
          Z = c.Y * Math.Sin(r) + c.Z * Math.Cos(r) }
    | None -> c

let rotateAxisY (c : Coords) (rad : float option) : Coords =
    match rad with
    | Some r ->
        { X = c.X * Math.Cos(r) + c.Z * Math.Sin(r)
          Y = c.Y
          Z = c.Z * Math.Cos(r) - c.X * Math.Sin(r) }
    | None -> c

let rotateAxisZ (c : Coords) (rad : float option) : Coords =
    match rad with
    | Some r ->
        { X = c.X * Math.Cos(r) - c.Y * Math.Sin(r)
          Y = c.X * Math.Sin(r) + c.Y * Math.Cos(r)
          Z = c.Z }
    | None -> c

// ============================================================================
// Write SVG.
// ============================================================================
let writeAtomStyle (sb : System.Text.StringBuilder) ((idx, _, _) : AtomInfo) =
    sb.Append($"\n.atom-{idx}{{fill:url(#radial-gradient-{idx});}}") |> ignore

let writeAtomDefs
    (sb : System.Text.StringBuilder)
    (pov : Coords)
    (unitSize : float)
    ((idx, atom, coords)
    : AtomInfo) =
    let distZ = abs (pov.Z - coords.Z)
    let ratio = unitSize / distZ
    let radius = atom.Radius * ratio
    let d1, d2, d3, d4, d5 = atom.Color.Gradient

    let opstr = Operators.string

    let cx = opstr coords.X
    let cy = opstr coords.Y
    let r = opstr (radius + 0.4)

    let dr1 = opstr dr1
    let dr2 = opstr dr2
    let dr3 = opstr dr3
    let dr4 = opstr dr4
    let dr5 = opstr dr5

    let d1r, d1g, d1b = opstr d1.R, opstr d1.G, opstr d1.B
    let d2r, d2g, d2b = opstr d2.R, opstr d2.G, opstr d2.B
    let d3r, d3g, d3b = opstr d3.R, opstr d3.G, opstr d3.B
    let d4r, d4g, d4b = opstr d4.R, opstr d4.G, opstr d4.B
    let d5r, d5g, d5b = opstr d5.R, opstr d5.G, opstr d5.B

    let atomDefsDescr =
        $"\n<radialGradient\
        \n\tid=\"radial-gradient-{idx}\"\
        \n\tcx=\"{cx}\"\
        \n\tcy=\"{cy}\"\
        \n\tfx=\"{cx}\"\
        \n\tfy=\"{cy}\"\
        \n\tr=\"{r}\"\
        \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{dr1}\" stop-color=\"rgb({d1r},{d1g},{d1b})\"/>\
        \n<stop offset=\"{dr2}\" stop-color=\"rgb({d2r},{d2g},{d2b})\"/>\
        \n<stop offset=\"{dr3}\" stop-color=\"rgb({d3r},{d3g},{d3b})\"/>\
        \n<stop offset=\"{dr4}\" stop-color=\"rgb({d4r},{d4g},{d4b})\"/>\
        \n<stop offset=\"{dr5}\" stop-color=\"rgb({d5r},{d5g},{d5b})\"/>\
        \n</radialGradient>"
    sb.Append(atomDefsDescr) |> ignore

let writeAtom
    (sb : System.Text.StringBuilder)
    (pov : Coords)
    (unitSize : float)
    ((idx, atom, coords)
    : AtomInfo) =
    let distZ = abs (pov.Z - coords.Z)
    let ratio = unitSize / distZ
    let radius = atom.Radius * ratio

    let opstr = Operators.string
    let cx = opstr coords.X
    let cy = opstr coords.Y
    let r = opstr radius

    let atomDescr =
        $"\n<circle\
        \n\tclass=\"atom-{idx}\"\
        \n\tcx=\"{cx}\"\
        \n\tcy=\"{cy}\"\
        \n\tr=\"{r}\"\
        \n/>"
    sb.Append(atomDescr) |> ignore

let writeSVG
    ((xMin, yMin, width, height) : ViewBox)
    (pov : Coords)
    (unitSize : float)
    (atoms : AtomInfo array)
    : string =
    let sb = System.Text.StringBuilder()
    let header =
        $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        \n<svg\
        \n\tid=\"Layer_1\"\
        \n\txmlns=\"http://www.w3.org/2000/svg\"\
        \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
        \n>"
    sb.Append(header) |> ignore
    sb.Append("\n<defs>\n<style>") |> ignore
    atoms |> Array.map (writeAtomStyle sb) |> ignore
    sb.Append("\n</style>") |> ignore
    atoms |> Array.map (writeAtomDefs sb pov unitSize) |> ignore
    sb.Append("\n</defs>") |> ignore
    atoms |> Array.map (writeAtom sb pov unitSize) |> ignore
    sb.Append("\n</svg>") |> ignore
    sb.ToString()

// ============================================================================
// Main.
// ============================================================================
let draw sdf =
    // Settings.
    let filterHydrogens : bool = false
    let numSteps : float = 20.0
    let rotation : (Coords -> float option -> Coords) = rotateAxisY

    // Parse atoms from SDF/Mol V2000 file; filter out hydrogen atoms.
    let filterAtoms (atomType : Atom) (atoms : AtomInfo array) : AtomInfo array =
        Array.filter (fun ((_, a, _) : AtomInfo) -> a <> atomType) atoms

    let atoms =
        parse_sdf sdf
        |> (if filterHydrogens then filterAtoms H else (fun arr -> arr))

    // Determine dimensions view box.
    let maxX = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.X) atoms |> Array.max
    let maxY = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.Y) atoms |> Array.max
    let maxZ = Array.map (fun ((_, _, c) : AtomInfo) -> abs c.Z) atoms |> Array.max
    let unitSize = (List.max [ maxX; maxY; maxZ ]) * 2.0
    let pov = { X = 0.0; Y = 0.0; Z = unitSize }
    let viewBox = (-unitSize, -unitSize, 2.0 * unitSize, 2.0 * unitSize)

//    // Generate SVGs.
//    for step in [ 1.0 .. 1.0 .. numSteps ] do
//        // Calculate rotation angle.
//        let rad : float option = Some ((step / numSteps) * 2.0 * Math.PI)
//
//        // Parse and prep atom data.
//        let rotatedAtoms =
//            atoms
//            |> Array.map (fun ((i, a, c) : AtomInfo) -> (i, a, rotation c rad))
//            |> Array.sortBy (fun ((_, _, c) : AtomInfo) -> - (abs (pov.Z - c.Z)))
//
//        // Write out atoms to SVG.
    let svg = writeSVG viewBox pov unitSize atoms
    svg



// ============================================================================
// Other.
// ============================================================================
module Storage =
    let todos = ResizeArray()

    let addTodo (todo: Todo) =
        if Todo.isValid todo.Description then
            todos.Add todo
            Ok()
        else
            Error "Invalid todo"

    do
        addTodo (Todo.create "Create new SAFE project")
        |> ignore

        addTodo (Todo.create "Write your app") |> ignore
        addTodo (Todo.create "Ship it !!!") |> ignore


let toBase64String (toEncode : string) : string =
    let bytes = System.Text.UTF8Encoding.GetEncoding(28591).GetBytes(toEncode)
    Convert.ToBase64String(bytes)

let todosApi =
    { getTodos = fun () -> async { return Storage.todos |> List.ofSeq }
      addTodo =
        fun todo ->
            async {
                return
                    match Storage.addTodo todo with
                    | Ok () -> todo
                    | Error e -> failwith e
            }
      render = fun sdf -> async {
          return sdf |> draw |> toBase64String
      }
    }

let webApp =
    Remoting.createApi ()
    |> Remoting.withRouteBuilder Route.builder
    |> Remoting.fromValue todosApi
    |> Remoting.buildHttpHandler

let app =
    application {
        url "http://*:8085"
        use_router webApp
        memory_cache
        use_static "public"
        use_gzip
    }

[<EntryPoint>]
let main _ =
    run app
    0