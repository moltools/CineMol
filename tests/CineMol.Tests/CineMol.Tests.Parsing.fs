module CineMol.Tests.Parsing

open CineMol.Types.Chem
open NUnit.Framework
open CineMol.Types.Fundamentals
open CineMol.Style
open CineMol.Parsing
    
[<TestFixture>]
type LinterTests () =
    
    [<Test>]
    member _.``Test parsing water molecule from correct string in SDF format.`` () =
        let src = """CT1000292221


  3  2  0  0  0               999 V2000
    0.0021   -0.0041    0.0020 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0110    0.9628    0.0073 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8669    1.3681    0.0011 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
        let expected: Molecule =
            { Atoms = [ Atom3D ({ Index = Index 1; Type = H; Color = CPK.Color H }, { X =  0.0021; Y = -0.0041; Z = 0.0020 }, PubChem.Radius H)
                        Atom3D ({ Index = Index 2; Type = O; Color = CPK.Color O }, { X = -0.0110; Y =  0.9628; Z = 0.0073 }, PubChem.Radius O)
                        Atom3D ({ Index = Index 3; Type = H; Color = CPK.Color H }, { X =  0.8669; Y =  1.3681; Z = 0.0011 }, PubChem.Radius H) ];
              Bonds = [ Bond { Index = Index 1; BeginAtomIndex = Index 1; EndAtomIndex = Index 2; Type = Single; Color = None }
                        Bond { Index = Index 2; BeginAtomIndex = Index 2; EndAtomIndex = Index 3; Type = Single; Color = None } ] }
        
        match (FileParser Sdf).Parse src with
        | Some [ result ] ->
            let atomsComparison = List.zip result.Atoms expected.Atoms |> List.map (fun (x, y) -> x = y) |> List.contains false
            let bondsComparison = List.zip result.Bonds expected.Bonds |> List.map (fun (x, y) -> x = y) |> List.contains false 
            match atomsComparison, bondsComparison with | false, false -> Assert.Pass() | _ -> Assert.Fail()
        | _ -> Assert.Fail()
       