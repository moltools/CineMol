module Server.Tests

open Expecto

open Shared
open Server
open Cinemol.Types

let server = testList "Server" [
    testCase "Calculate dot product of two vectors" <| fun _ ->
        let u: Vector = { X = 1.0; Y = 2.0; Z = 3.0 }
        let v: Vector = { X = 4.0; Y = -5.0; Z = 6.0 }
        let result = u.Dot v
        let expected = 12.0
        Expect.equal expected result $"{expected} is not {result}"

    testCase "Calculate cross product of two vectors" <| fun _ ->
        let u: Vector = { X = 3.0; Y = -3.0; Z = 1.0 }
        let v: Vector = { X = 4.0; Y = 9.0; Z = 2.0 }
        let result = u.Cross v
        let expected: Vector = { X = -15.0; Y = -2.0; Z = 39.0 }
        Expect.equal expected result $"{expected} is not {result}"
]

let all =
    testList "All"
        [
            Shared.Tests.shared
            server
        ]

[<EntryPoint>]
let main _ = runTestsWithCLIArgs [] [||] all