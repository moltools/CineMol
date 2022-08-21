// Author: David Meijer
// Description: 3D SVG drawer for small molecules.
// Usage: dotnet fsi cinemol.fsx <path to SDF/Mol V2000 file>
open System
open System.IO
open System.Text.RegularExpressions

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
let args : string array = fsi.CommandLineArgs |> Array.tail
let input_path : string = args.[0]
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

let parse_sdf (sdf_path : string) : AtomInfo array =
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

let rotateAxisY (c : Coords) (rad : float option) : Coords = 
    match rad with 
    | Some r ->
        { X = c.X * Math.Cos(r) + c.Z * Math.Sin(r)
          Y = c.Y
          Z = c.Z * Math.Cos(r) - c.X * Math.Sin(r) }   
    | None -> c 

// ============================================================================
// Write SVG.
// ============================================================================
let writeAtomStyle (w : StreamWriter) ((idx, atom, _) : AtomInfo) =
    w.WriteLine($".atom-{idx}{{fill:url(#radial-gradient-{idx});}}")

let writeAtomDefs 
    (w : StreamWriter) 
    (pov : Coords)
    (unitSize : float)
    ((idx, atom, coords) 
    : AtomInfo) =
    let distZ = abs (pov.Z - coords.Z)
    let ratio = unitSize / distZ
    let radius = atom.Radius * ratio
    let d1, d2, d3, d4, d5 = atom.Color.Gradient
    let atomDefsDescr = 
        $"<radialGradient\
        \n\tid=\"radial-gradient-{idx}\"\
        \n\tcx=\"{coords.X}\"\
        \n\tcy=\"{coords.Y}\"\
        \n\tfx=\"{coords.X}\"\
        \n\tfy=\"{coords.Y}\"\
        \n\tr=\"{radius + 0.4}\"\
        \n\tgradientTransform=\"matrix(1, 0, 0, 1, 0, 0)\"\
        \n\tgradientUnits=\"userSpaceOnUse\"\
        \n>\
        \n<stop offset=\"{dr1}\" stop-color=\"rgb({d1.R},{d1.G},{d1.B})\"/>\
        \n<stop offset=\"{dr2}\" stop-color=\"rgb({d2.R},{d2.G},{d2.B})\"/>\
        \n<stop offset=\"{dr3}\" stop-color=\"rgb({d3.R},{d3.G},{d3.B})\"/>\
        \n<stop offset=\"{dr4}\" stop-color=\"rgb({d4.R},{d4.G},{d4.B})\"/>\
        \n<stop offset=\"{dr5}\" stop-color=\"rgb({d5.R},{d5.G},{d5.B})\"/>\
        \n</radialGradient>"
    w.WriteLine(atomDefsDescr)

let writeAtom 
    (w : StreamWriter) 
    (pov : Coords) 
    (unitSize : float) 
    ((idx, atom, coords) 
    : AtomInfo) =
    let distZ = abs (pov.Z - coords.Z)
    let ratio = unitSize / distZ
    let radius = atom.Radius * ratio
    let atomDescr = 
        $"<circle\
        \n\tclass=\"atom-{idx}\"\
        \n\tcx=\"{coords.X}\"\
        \n\tcy=\"{coords.Y}\"\
        \n\tr=\"{radius}\"\
        \n/>"
    w.WriteLine(atomDescr)

let writeSVG 
    (out_path : string)
    ((xMin, yMin, width, height) : ViewBox) 
    (pov : Coords)
    (unitSize : float)
    (atoms : AtomInfo array) =
    use w = new StreamWriter(out_path)
    let header =
        $"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
        \n<svg\
        \n\tid=\"Layer_1\"\
        \n\txmlns=\"http://www.w3.org/2000/svg\"\
        \n\tviewBox=\"{xMin} {yMin} {width} {height}\"\
        \n>"
    w.WriteLine(header)
    w.WriteLine("<defs>\n<style>")
    atoms |> Array.map (writeAtomStyle w) |> ignore
    w.WriteLine("</style>")
    atoms |> Array.map (writeAtomDefs w pov unitSize) |> ignore
    w.WriteLine("</defs>")
    atoms |> Array.map (writeAtom w pov unitSize) |> ignore 
    w.WriteLine("</svg>")

// ============================================================================
// Main.
// ============================================================================
let main =
    // Settings.
    let filterHydrogens : bool = false
    let unitSize = 10.0
    let pov = { X = 0.0; Y = 0.0; Z = unitSize }
    let viewBox = (-unitSize, -unitSize, 2.0 * unitSize, 2.0 * unitSize)

    let numSteps = 20.0
    for step in [ 1.0 .. 1.0 .. numSteps ] do
        let rad : float option = Some ((step / numSteps) * 2.0 * Math.PI)

        // Parse and prep atom data.
        let filterAtoms (atomType : Atom) (atoms : AtomInfo array) : AtomInfo array =
            Array.filter (fun ((_, a, _) : AtomInfo) -> a <> atomType) atoms
        let atoms = 
            parse_sdf input_path
            |> (if filterHydrogens then filterAtoms H else (fun arr -> arr))
            |> Array.map (fun ((i, a, c) : AtomInfo) -> (i, a, rotateAxisY c rad))
            |> Array.sortBy (fun ((_, _, c) : AtomInfo) -> - (abs (pov.Z - c.Z)))
        
        // Write out atoms to SVG.
        writeSVG $"out/mol_{step}.svg" viewBox pov unitSize atoms |> ignore

    // Exit code.
    0

[<EntryPoint>]
main |> printfn "%A"
