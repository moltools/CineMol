namespace YourNamespace

module YourModule =

    // Define a function that you want to call from Python
    let add x y =
        x + y

    // Define another function
    let subtract x y =
        x - y

    // Define a type with methods (you can access methods from Python too)
    type Calculator() =
        member this.add x y =
            x + y
        member this.subtract x y =
            x - y
