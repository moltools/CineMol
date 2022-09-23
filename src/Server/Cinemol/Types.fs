namespace Cinemol

module Types =

    let dr1, dr2, dr3, dr4, dr5 = 0.0, 0.11, 0.34, 0.66, 1.0  // Diffusion rates

    type Gradient = Color * Color * Color * Color * Color

    and Color =
        { R : int
          G : int
          B : int }
        member this.Diffuse (factor : float) : Color =
            { R = int ((float this.R) * factor)
              G = int ((float this.G) * factor)
              B = int ((float this.B) * factor) }
        member this.Gradient : Gradient =
            ( this.Diffuse (1.0 - dr1),
              this.Diffuse (1.0 - dr2),
              this.Diffuse (1.0 - dr3),
              this.Diffuse (1.0 - dr4),
              this.Diffuse (1.0 - dr5) )

    type Index = int

    type Atom =
        | C | N | O | S | H
        member this.Radius : float =
            match this with
            | H -> 1.0
            | _ -> 1.2
        member this.Color : Color =
            match this with
            | C -> { R = 128; G = 128; B = 128 }
            | N -> { R = 0; G = 0; B = 255 }
            | O -> { R = 255; G = 0; B = 0 }
            | S -> { R = 255; G = 255; B = 0 }
            | H -> { R = 220; G = 220; B = 220 }

    type Coords =
        { X : float
          Y : float
          Z : float }
        static member (-) (c1 : Coords, c2 : Coords) : Coords =
            { X = c1.X - c2.X; Y = c1.Y - c2.Y; Z = c1.Z - c2.Z }
        static member Pow (c : Coords) (d : float) : Coords =
            { X = c.X ** d; Y = c.Y ** d; Z = c.Z ** d }
        static member Sum (c : Coords) : float =
            c.X + c.Y + c.Z

    type AtomInfo = Index * Atom * Coords

    type ViewBox = float * float * float * float

    type Depiction =
        | Filled
        | BallAndStick
