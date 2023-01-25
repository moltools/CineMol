module CineMol.Helpers

open System

/// Round floating point number to desired number of digits.
let round (digits: int) (v: float) = Math.Round(v, digits)

/// Take absolute value of floating point number.
let abs (v: float) = Math.Sqrt(v ** 2.0)

/// Clamp a floating point value between two other floating point values.
let clamp lowerBound upperBound (v: float) =
    if v < lowerBound then lowerBound
    elif v > upperBound then upperBound
    else v
    
/// Try to cast string.
let inline tryCast (cast: string -> 'a) s =
    try s |> cast |> Some 
    with :? FormatException -> None