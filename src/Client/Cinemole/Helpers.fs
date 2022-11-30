module Client.Cinemole.Helpers

open System
open Encoder


let floatToStr (f: float) : string = Operators.string f
let intToStr (d: int) : string = Operators.string d

let round (digits: int) (f: float) = Math.Round(f, digits)

let abs f = (f ** 2.0) ** 0.5

let toBase64String (toEncode : string) : string =
//    let bytes = System.Text.UTF8Encoding.GetEncoding(28591).GetBytes(toEncode)
    let bytes = encodeString toEncode
    Convert.ToBase64String(bytes)