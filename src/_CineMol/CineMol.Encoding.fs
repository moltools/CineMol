module CineMol.Encoding

type Encoding = | ISO_8859_1
    with
    member this.EncodeChar c =
        // Pick encoding type.
        match this with

        // Base64 encode character according to ISO 8859-1.
        | ISO_8859_1 ->
            match c with
            | '\t'  | '\n'  | ' '  -> Some 32uy
            | '!'  -> Some 33uy  | '"'  -> Some 34uy  | '#'  -> Some 35uy
            | '$'  -> Some 36uy  | '%'  -> Some 37uy  | '&'  -> Some 38uy
            | '''  -> Some 39uy  | '('  -> Some 40uy  | ')'  -> Some 41uy
            | '*'  -> Some 42uy  | '+'  -> Some 43uy  | ','  -> Some 44uy
            | '-'  -> Some 45uy  | '.'  -> Some 46uy  | '/'  -> Some 47uy
            | '0'  -> Some 48uy  | '1'  -> Some 49uy  | '2'  -> Some 50uy
            | '3'  -> Some 51uy  | '4'  -> Some 52uy  | '5'  -> Some 53uy
            | '6'  -> Some 54uy  | '7'  -> Some 55uy  | '8'  -> Some 56uy
            | '9'  -> Some 57uy  | ':'  -> Some 58uy  | ';'  -> Some 59uy
            | '<'  -> Some 60uy  | '='  -> Some 61uy  | '>'  -> Some 62uy
            | '?'  -> Some 63uy  | '@'  -> Some 64uy  | 'A'  -> Some 65uy
            | 'B'  -> Some 66uy  | 'C'  -> Some 67uy  | 'D'  -> Some 68uy
            | 'E'  -> Some 69uy  | 'F'  -> Some 70uy  | 'G'  -> Some 71uy
            | 'H'  -> Some 72uy  | 'I'  -> Some 73uy  | 'J'  -> Some 74uy
            | 'K'  -> Some 75uy  | 'L'  -> Some 76uy  | 'M'  -> Some 77uy
            | 'N'  -> Some 78uy  | 'O'  -> Some 79uy  | 'P'  -> Some 80uy
            | 'Q'  -> Some 81uy  | 'R'  -> Some 82uy  | 'S'  -> Some 83uy
            | 'T'  -> Some 84uy  | 'U'  -> Some 85uy  | 'V'  -> Some 86uy
            | 'W'  -> Some 87uy  | 'X'  -> Some 88uy  | 'Y'  -> Some 89uy
            | 'Z'  -> Some 90uy  | '['  -> Some 91uy  | '\\' -> Some 92uy
            | ']'  -> Some 93uy  | '^'  -> Some 94uy  | '_'  -> Some 95uy
            | '`'  -> Some 96uy  | 'a'  -> Some 97uy  | 'b'  -> Some 98uy
            | 'c'  -> Some 99uy  | 'd'  -> Some 100uy | 'e'  -> Some 101uy
            | 'f'  -> Some 102uy | 'g'  -> Some 103uy | 'h'  -> Some 104uy
            | 'i'  -> Some 105uy | 'j'  -> Some 106uy | 'k'  -> Some 107uy
            | 'l'  -> Some 108uy | 'm'  -> Some 109uy | 'n'  -> Some 110uy
            | 'o'  -> Some 111uy | 'p'  -> Some 112uy | 'q'  -> Some 113uy
            | 'r'  -> Some 114uy | 's'  -> Some 115uy | 't'  -> Some 116uy
            | 'u'  -> Some 117uy | 'v'  -> Some 118uy | 'w'  -> Some 119uy
            | 'x'  -> Some 120uy | 'y'  -> Some 121uy | 'z'  -> Some 122uy
            | '{'  -> Some 123uy | '|'  -> Some 124uy | '}'  -> Some 125uy
            | '~'  -> Some 126uy
            | c ->
                printf $"Unexpected char for ISO 8859-1 encoding: '{c}'"
                None

    member this.Encode (src: string) =
        src
        |> Seq.toArray
        |> Array.map (fun c -> this.EncodeChar c)
        |> Array.choose id // We drop unexpected characters.