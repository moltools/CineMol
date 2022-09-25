module Server.Tests

open Expecto

open Shared
open Server
open Cinemol.Helpers
open Cinemol.Types
open Cinemol.Geometry

let server = testList "Server" [
    testCase "Calculate dot product of two vectors" <| fun _ ->
        let u: Vector = { X = 1.0; Y = 2.0; Z = 3.0 }
        let v: Vector = { X = 4.0; Y = -5.0; Z = 6.0 }
        Expect.equal 12.0 (u.Dot v) "failed"

    testCase "Calculate cross product of two vectors" <| fun _ ->
        let u: Vector = { X = 3.0; Y = -3.0; Z = 1.0 }
        let v: Vector = { X = 4.0; Y = 9.0; Z = 2.0 }
        let expected: Vector = { X = -15.0; Y = -2.0; Z = 39.0 }
        Expect.equal expected (u.Cross v) "failed"

    testCase "Calculate plane of intersection of two atoms with the same center and radius" <| fun _ ->
        let a1 = createAtom 1 C { X = 0.0; Y = 0.0; Z = 0.0 } 1.0
        let a2 = createAtom 2 C { X = 0.0; Y = 0.0; Z = 0.0 } 1.0
        Expect.equal NoIntersection (a1.Intersect a2) "failed"

    testCase "Calculate plane of intersection of two atoms that have intersection point" <| fun _ ->
        let a1 = createAtom 1 C { X = 0.0; Y = 0.0; Z = 0.0 } 1.0
        let a2 = createAtom 2 C { X = 2.0; Y = 0.0; Z = 0.0 } 1.0
        let result = a1.Intersect a2
        Expect.equal (IntersectionPoint { X = 1.0; Y = 0.0; Z = 0.0 }) result "failed"

    testCase "Calculate plane of intersection of two atoms that have no intersection" <| fun _ ->
        let a1 = createAtom 1 C { X = 0.0; Y = 0.0; Z = 0.0 } 0.5
        let a2 = createAtom 2 C { X = 2.0; Y = 0.0; Z = 0.0 } 0.5
        let result = a1.Intersect a2
        Expect.equal NoIntersection result "failed"

    testCase "Calculate plane of intersection of two atoms that have intersection circle" <| fun _ ->
        let a1 = createAtom 1 C { X = 0.0; Y = 0.0; Z = 0.0 } 1.5
        let a2 = createAtom 2 C { X = 2.0; Y = 0.0; Z = 0.0 } 1.5
        let result = a1.Intersect a2
        let cExpected: Point = { X = 1.0; Y = 0.0; Z = 0.0 }
        let rExpected: Radius = 1.118
        let vExpected: Vector = { X = 2.0; Y = 0.0; Z = 0.0 }
        match result with
        | IntersectionCircle (c, r, v) ->
            Expect.equal cExpected c "failed"
            Expect.equal rExpected (round 3 r) "failed"
            Expect.equal vExpected v "failed"
        | _ -> Expect.isTrue false "failed"

    testCase "Calculate intersection points of two circles" <| fun _ ->
        let c1, r1: Point * Radius = { X = 0.0; Y = 0.0; Z = 0.0 }, 2.0
        let c2, r2: Point * Radius = { X = 2.0; Y = 0.0; Z = 0.0 }, 2.0
        match intersectionCircles c1 r1 c2 r2 with
        | None -> Expect.isTrue false "failed"
        | Some (p1, p2) ->
            Expect.equal 1.0 p1.X "failed"
            Expect.equal -1.732 (round 3 p1.Y) "failed"
            Expect.equal 0.0 p1.Z "failed"
            Expect.equal 1.0 p2.X "failed"
            Expect.equal 1.732 (round 3 p2.Y) "failed"
            Expect.equal 0.0 p2.Z "failed"
]

let all = testList "All" [ Shared.Tests.shared; server ]

[<EntryPoint>]
let main _ = runTestsWithCLIArgs [] [||] all