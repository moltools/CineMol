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

        // let atomIndicesToHighlight = [ 2; 6; 8; 11; 12 ]
        // let highlightColor = Types.Color (238, 130, 238) // Violet 
        // let defaultColor = Types.Color (128, 128, 128) // Grey

        // let molecule = 
        //     { molecule with 
        //         Atoms = 
        //             molecule.Atoms 
        //             |> List.map (fun atom ->
        //                 if List.contains atom.Index atomIndicesToHighlight then
        //                     { atom with Color = highlightColor }
        //                 else
        //                     { atom with Color = defaultColor }) }

        // let rotationAxis = Types.Geometry.Axis.X
        // let rotationRads = Math.PI / 2.0 // 90 degrees
        // let molecule = Drawing.rotate molecule rotationAxis rotationRads

        let options : Types.Drawing.DrawingOptions = 
            { ViewBox = None 
              Style = Types.ModelStyle.SpaceFilling  
              DisplayHydrogenAtoms = false
              Resolution = 200 }
              
        let svg = Drawing.draw molecule options |> fst
        System.IO.File.WriteAllText(outFilePath, svg.ToString())

        0

    | None -> 1

[<EntryPoint>]
Environment.GetCommandLineArgs() |> main |> printfn "%A"