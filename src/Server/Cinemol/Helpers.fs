namespace Cinemol

module Helpers =

    open System

    let toBase64String (toEncode : string) : string =
        let bytes = System.Text.UTF8Encoding.GetEncoding(28591).GetBytes(toEncode)
        Convert.ToBase64String(bytes)