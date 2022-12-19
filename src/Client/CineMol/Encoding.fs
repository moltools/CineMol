module Client.CineMol.Encoding

open System

exception UnexpectedCharToEncode of string

type Encoding = | ISO_8859_1
    with
    member x.GetEncoding (c: char) : byte =
        match x with

        // Encode character according to ISO 8859-1.
        | ISO_8859_1 ->
            match c with
            | '\t'  | '\n'  | ' '  -> 32uy
            | '!'  -> 33uy  | '"'  -> 34uy  | '#'  -> 35uy  | '$'  -> 36uy
            | '%'  -> 37uy  | '&'  -> 38uy  | '''  -> 39uy  | '('  -> 40uy
            | ')'  -> 41uy  | '*'  -> 42uy  | '+'  -> 43uy  | ','  -> 44uy
            | '-'  -> 45uy  | '.'  -> 46uy  | '/'  -> 47uy  | '0'  -> 48uy
            | '1'  -> 49uy  | '2'  -> 50uy  | '3'  -> 51uy  | '4'  -> 52uy
            | '5'  -> 53uy  | '6'  -> 54uy  | '7'  -> 55uy  | '8'  -> 56uy
            | '9'  -> 57uy  | ':'  -> 58uy  | ';'  -> 59uy  | '<'  -> 60uy
            | '='  -> 61uy  | '>'  -> 62uy  | '?'  -> 63uy  | '@'  -> 64uy
            | 'A'  -> 65uy  | 'B'  -> 66uy  | 'C'  -> 67uy  | 'D'  -> 68uy
            | 'E'  -> 69uy  | 'F'  -> 70uy  | 'G'  -> 71uy  | 'H'  -> 72uy
            | 'I'  -> 73uy  | 'J'  -> 74uy  | 'K'  -> 75uy  | 'L'  -> 76uy
            | 'M'  -> 77uy  | 'N'  -> 78uy  | 'O'  -> 79uy  | 'P'  -> 80uy
            | 'Q'  -> 81uy  | 'R'  -> 82uy  | 'S'  -> 83uy  | 'T'  -> 84uy
            | 'U'  -> 85uy  | 'V'  -> 86uy  | 'W'  -> 87uy  | 'X'  -> 88uy
            | 'Y'  -> 89uy  | 'Z'  -> 90uy  | '['  -> 91uy  | '\\' -> 92uy
            | ']'  -> 93uy  | '^'  -> 94uy  | '_'  -> 95uy  | '`'  -> 96uy
            | 'a'  -> 97uy  | 'b'  -> 98uy  | 'c'  -> 99uy  | 'd'  -> 100uy
            | 'e' -> 101uy  | 'f'  -> 102uy | 'g'  -> 103uy | 'h'  -> 104uy
            | 'i' -> 105uy  | 'j'  -> 106uy | 'k'  -> 107uy | 'l'  -> 108uy
            | 'm' -> 109uy  | 'n'  -> 110uy | 'o'  -> 111uy | 'p'  -> 112uy
            | 'q' -> 113uy  | 'r'  -> 114uy | 's'  -> 115uy | 't'  -> 116uy
            | 'u' -> 117uy  | 'v'  -> 118uy | 'w'  -> 119uy | 'x'  -> 120uy
            | 'y' -> 121uy  | 'z'  -> 122uy | '{'  -> 123uy | '|'  -> 124uy
            | '}' -> 125uy  | '~'  -> 126uy
            | _   ->
                $"Unknown character to encode {x}: {c}"
                |> UnexpectedCharToEncode
                |> raise

let encode (encoding: Encoding) (src: string) : byte[] =
    src
    |> Seq.toArray
    |> Array.map (fun c -> encoding.GetEncoding c)

let toBase64String (toEncode : string) : string =
    let bytes =  encode ISO_8859_1 toEncode
    Convert.ToBase64String(bytes)
