module Client.CineMol.Helpers

open System

let floatToStr (f: float) : string = Operators.string f
let intToStr (d: int) : string = Operators.string d

let round (n: int) (f: float) = Math.Round(f, n)

let abs f = (f ** 2.0) ** 0.5