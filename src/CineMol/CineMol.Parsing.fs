module CineMol.Parsing

open System.Text.RegularExpressions
open Helpers
open Style 
open Types.Geometry
open Types.Chem

/// <summary>
/// Supported file types for parsing molecule structures from.
/// </summary>
exception ParserError of string 

/// <summary>
/// Supported file types for parsing molecule structures from.
/// </summary>
type FileType = | Sdf

/// <summary>
/// Parse input file containing molecule structures.
/// </summary> 
type FileParser = FileParser of FileType 
    with
    member this.Parse (fileContent: string) : Molecule list option =
        match (match this with | FileParser fileType -> fileType) with
        
        // TODO: check if file type corresponds with file extension and contents of file
        // TODO: limit file size for uploading
        
        // Parse molecules from SDF file.
        | Sdf ->
            /// Defines regex for whitespace.
            let s = @"\s{1,}"
            
            /// Defines regex for floating point value.
            let d = @"[-+]?[0-9]*\.?[0-9]+"
            
            /// Defines regex for captured floating point value.
            let d_cap = $"({d})"
            
            /// Defines regex for captured word.
            let w_cap = @"(\w+)"
            
            /// Constructs regex for atom line in input file.
            let atomLine : string =
                [ "$" ]
                |> (@) [ for _ in [ 0 .. 11 ] do yield s + d ]
                |> (@) [ s + w_cap]
                |> (@) [ for _ in [ 0 .. 2 ] do yield s + d_cap ]
                |> (@) [ "^" ]
                |> String.concat ""
                
            /// Constructs regex for bond line in input file.
            let bondLine : string =
                [ "$" ]
                |> (@) [ for _ in [ 0 .. 3 ] do yield s + d ]
                |> (@) [ for _ in [ 0 .. 2 ] do yield s + d_cap ]
                |> (@) [ "^" ]
                |> String.concat ""
            
            /// Active pattern for matching and parsing atom line in input file.
            let (|AtomLine|_|) (line: string) : string list option =
                let m = Regex.Match(line, atomLine)
                if m.Success then Some(List.tail [ for g in m.Groups -> g.Value ])
                else None
            
            /// Constructs regex for matching and parsing bond line in input file.
            let (|BondLine|_|) (line: string) : string list option =
                let m = Regex.Match(line, bondLine)
                if m.Success then Some(List.tail [ for g in m.Groups -> g.Value ])
                else None
            
            /// Parse input file.
            let mutable count = 0
            let mutable atoms: Atom list = []
            let mutable bonds: Bond list = []
            
            [
                for line in fileContent.Split [| '\n' |] do
                    match line with
                    
                    // Molecule structure delimiter.
                    | line when line.Contains("$$$$") = true ->
                        yield { Atoms = atoms; Bonds = bonds }
                        
                        count <- count + 1
                        atoms <- []
                        bonds <- []
                    
                    // Atom line.
                    | AtomLine [ x; y; z; atomSymbol ] ->
                        match
                            tryCast float x,
                            tryCast float y,
                            tryCast float z,
                            AtomType.FromString atomSymbol
                            with
                        
                        // Able to cast all data to appropriate types.
                        | Some x, Some y, Some z, Some atomType ->
                            let atom : Atom =
                                { Index = atoms.Length + (bonds.Length * 2) + 1
                                  Type = atomType
                                  Color = CPK.Color atomType
                                  Opacity = 1.0
                                  Position = { X = x; Y = y; Z = z }
                                  Radius = PubChem.Radius atomType }

                            atoms <- atoms @ [ atom ]
                        
                        // Unable to cast all data to appropriate types.
                        | _ ->
                            $"Unable to parse atom line: '{line}'"
                            |> ParserError |> raise 
                    
                    // Bond line.
                    | BondLine [ s_idx; e_idx; bondType ] ->
                        match
                            tryCast int s_idx,
                            tryCast int e_idx,
                            BondType.FromString bondType
                            with
                        
                        // Able to cast all data to appropriate types.
                        | Some s_idx, Some e_idx, Some bondType ->
                            let bond : Bond =
                                { BeginIndex = atoms.Length + (bonds.Length * 2) + 1
                                  EndIndex = atoms.Length + (bonds.Length * 2) + 2
                                  Type = bondType
                                  BeginAtomIndex = s_idx
                                  EndAtomIndex = e_idx
                                  Opacity = None 
                                  Color = None
                                  Radius = 0.5 }

                            bonds <- bonds @ [ bond ]  
                        
                        // Unable to cast all data to appropriate types.
                        | _ ->
                            $"Unable to parse bond line: '{line}'"
                            |> ParserError |> raise 
                    
                    // Skip other lines.
                    | _ -> ()
                    
            ]
            |> List.map (fun mol -> mol.AdjustForCentroid())
            |> Some 