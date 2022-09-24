module Shared.Tests

#if FABLE_COMPILER
open Fable.Mocha
#else
open Expecto
#endif

open Shared

let shared = testList "Shared" [
    testCase "Dummy test Shared" <| fun _ ->
        Expect.equal true true "Should be true"
]