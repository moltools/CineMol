#r "../src/CineMol/bin/Debug/net6.0/CineMol.dll"

open System
open CineMol

let main (argv : string array) =
    let sdfFilePath = argv[2]
    let outFilePath = argv[3]

    let sdfFile = System.IO.File.ReadAllText(sdfFilePath)
    let parser = Parsing.FileParser Parsing.Sdf 

    match parser.Parse sdfFile with
    | Some molecules -> 
        let molecule = List.head molecules

        let rotationAxis = Types.Geometry.Axis.X
        let rotationRads = Math.PI / 2.0 // 90 degrees
        let molecule = Drawing.rotate molecule rotationAxis rotationRads

        let options : Types.Drawing.DrawingOptions = 
            { ViewBox = None 
              Style = Types.ModelStyle.SpaceFilling  
              DisplayHydrogenAtoms = true
              Resolution = 200 }
              
        let svg = Drawing.draw molecule options |> fst
        System.IO.File.WriteAllText(outFilePath, svg.ToString())

        0

    | None -> 1

[<EntryPoint>]
Environment.GetCommandLineArgs() |> main |> printfn "%A"