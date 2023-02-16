module CineMol.Helpers

open System

/// <summary>
/// Round floating point number to desired number of digits.
/// </summary>
let round (digits: int) (v: float) = Math.Round(v, digits)

/// <summary>
/// Take absolute value of floating point number.
/// </summary>
let abs (v: float) = Math.Sqrt(v ** 2.0)

/// <summary>
/// Clamp a floating point value between two other floating point values.
/// </summary>
let clamp lowerBound upperBound (v: float) =
    if v < lowerBound then lowerBound
    elif v > upperBound then upperBound
    else v
    
/// <summary>
/// Try to cast string.
/// </summary>
let inline tryCast (cast: string -> 'a) s =
    try s |> cast |> Some 
    with :? FormatException -> None
    
/// <summary>
/// Enumerate a list.
/// </summary>
let inline enumerate (items: 'a list) = List.zip [0 .. 1 .. items.Length] items 