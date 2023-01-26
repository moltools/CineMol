module CineMol.Tests.Parsing

open NUnit.Framework
open CineMol.Parsing
    
[<TestFixture>]
type CineMolTests () =
    
    [<Test>]
    member _.``Test parsing water molecule from correct string in SDF format`` () =
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
        match (FileParser Sdf).Parse src with
        | Some [ { Atoms = atoms; Bonds = bonds } ] ->
            match atoms, bonds with
            | atoms, bonds when atoms.Length = 3 && bonds.Length = 2 -> Assert.Pass()
            | _ -> Assert.Fail()
        | _ -> Assert.Fail()