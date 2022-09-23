namespace Cinemol

module Helpers =

    open System

    let round (digits: int) (f: float) = Math.Round(f, digits)
    let abs f = (f ** 2.0) ** 0.5

    let toBase64String (toEncode : string) : string =
        let bytes = System.Text.UTF8Encoding.GetEncoding(28591).GetBytes(toEncode)
        Convert.ToBase64String(bytes)