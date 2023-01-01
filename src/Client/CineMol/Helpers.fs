module Client.CineMol.Helpers

open System

let round (digits: int) (f: float) = Math.Round(f, digits)
let floatToStr (f: float) : string = Operators.string f
let intToStr (d: int) : string = Operators.string d
let abs f = (f ** 2.0) ** 0.5