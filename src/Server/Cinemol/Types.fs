namespace Cinemol

module Types =

    open System

    open Cinemol.Helpers

    type Index = int

    type Radius = float

    type DiffusionRate = float

    let diffusionRate1: DiffusionRate = 0.0
    let diffusionRate2: DiffusionRate = 0.11
    let diffusionRate3: DiffusionRate = 0.34
    let diffusionRate4: DiffusionRate = 0.66
    let diffusionRate5: DiffusionRate = 01.0

    type Gradient = Color * Color * Color * Color * Color

    and Color = { R: int; G: int; B: int }
        with
        member x.Diffuse (factor : float) : Color =
            { R = int ((float x.R) * factor)
              G = int ((float x.G) * factor)
              B = int ((float x.B) * factor) }

        member x.Gradient : Gradient =
            ( x.Diffuse (1.0 - diffusionRate1),
              x.Diffuse (1.0 - diffusionRate2),
              x.Diffuse (1.0 - diffusionRate3),
              x.Diffuse (1.0 - diffusionRate4),
              x.Diffuse (1.0 - diffusionRate5) )

    let dodgerBlue: Color = { R = 1; G = 122; B = 255 }
    let mutedRed: Color = { R = 215; G = 80; B = 77 }
    let grey: Color = { R = 104; G = 104; B = 104 }
    let lightGrey: Color = { R = 238 ; G = 238 ; B = 238 }
    let neon: Color = { R = 108; G = 71; B = 255 }
    let turquoise: Color = { R = 57; G = 192; B = 200 }
    let safetyOrange: Color = { R = 249; G = 99; B = 0 }
    let wildMelon: Color = { R = 243; G = 37; B = 113 }
    let coral: Color = { R = 255; G = 147; B = 130 }
    let amber: Color = { R = 245; G = 201; B = 0 }
    let tan: Color = { R = 205; G = 173; B = 122 }
    let celery: Color = { R = 170; G = 187; B = 93 }

    type Axis = | X | Y | Z
        with
        member x.RotationMatrix : Point * float -> Point =
            match x with
            | X ->
                (fun (p: Point, rad: float) ->
                    { X = p.X
                      Y = p.Y * Math.Cos(rad) - p.Z * Math.Sin(rad)
                      Z = p.Y * Math.Sin(rad) + p.Z * Math.Cos(rad) })
            | Y ->
                (fun (p: Point, rad: float) ->
                    { X = p.X * Math.Cos(rad) + p.Z * Math.Sin(rad)
                      Y = p.Y
                      Z = p.Z * Math.Cos(rad) - p.X * Math.Sin(rad) })
            | Z ->
                (fun (p: Point, rad: float) ->
                    { X = p.X * Math.Cos(rad) - p.Y * Math.Sin(rad)
                      Y = p.X * Math.Sin(rad) + p.Y * Math.Cos(rad)
                      Z = p.Z })

    and Point = { X: float; Y: float; Z: float }
        with
        static member (-) (p1: Point, c2: Point) : Point = { X = p1.X - c2.X; Y = p1.Y - c2.Y; Z = p1.Z - c2.Z }

        static member Pow (p: Point) (d: float) : Point = { X = p.X ** d; Y = p.Y ** d; Z = p.Z ** d }

        static member Sum (p: Point) : float = p.X + p.Y + p.Z

        member p1.Distance (p2: Point) : float = Math.Sqrt(Point.Sum(Point.Pow (p1 - p2) 2.0))

        member p.Rotate (axis: Axis) (rad: float) : Point = axis.RotationMatrix(p, rad)

        member p1.FindVector (p2: Point) : Vector = { X = p2.X - p1.X; Y = p2.Y - p1.Y; Z = p2.Z - p1.Z }

    and Vector = { X: float; Y: float; Z: float }
        with
        member u.SumOfSquares : float = u.X ** 2.0 + u.Y ** 2.0 + u.Z ** 2.0

        member u.Magnitude : float = Math.Sqrt(u.SumOfSquares)

        member u.Dot (v: Vector) : float = u.X * v.X + u.Y * v.Y + u.Z * v.Z

        member u.Cross (v: Vector) : Vector =
            { X = u.Y * v.Z - u.Z * v.Y
              Y = u.Z * v.X - u.X * v.Z
              Z = u.X * v.Y - u.Y * v.X }

        member x.ProjectVector (v: Vector) : float = (v.Dot x) / v.Magnitude

    type SphereSphereIntersection =
        | Eclipsed
        | NoIntersection
        | IntersectionPoint of Point
        | IntersectionCircle of Point * Radius * Vector

    type Atom = | C | N | O | S | H
        with
        member x.Radius : float =
            match x with
            | S -> 1.2
            | H -> 0.6
            | _ -> 1.0

        member x.Color : Color =
            match x with
            | C -> grey
            | N -> dodgerBlue
            | O -> mutedRed
            | S -> amber
            | H -> lightGrey

    type AtomInfo =
        { Index: Index
          Type: Atom
          OriginalCenter: Point
          ProjectedCenter: Point
          OriginalRadius: Radius
          ProjectedRadius: Radius }
        with
        member x.Rotate (axis: Axis) (rad: float) : AtomInfo = { x with OriginalCenter = x.OriginalCenter.Rotate axis rad }

        member this.Intersect (other: AtomInfo) : SphereSphereIntersection =
            let dist = this.OriginalCenter.Distance other.OriginalCenter

            match dist with
            | d when d > (this.OriginalRadius + other.OriginalRadius) ||
                     (d = 0.0 && this.OriginalRadius = other.OriginalRadius)
                      -> NoIntersection
            | d when (d + this.OriginalRadius) < other.OriginalRadius  -> Eclipsed
            | _ ->
                // Intersection plane
                let A = 2.0 * (other.OriginalCenter.X - this.OriginalCenter.X)
                let B = 2.0 * (other.OriginalCenter.Y - this.OriginalCenter.Y)
                let C = 2.0 * (other.OriginalCenter.Z - this.OriginalCenter.Z)
                let D = this.OriginalCenter.X ** 2.0 - other.OriginalCenter.X ** 2.0 +
                        this.OriginalCenter.Y ** 2.0 - other.OriginalCenter.Y ** 2.0 +
                        this.OriginalCenter.Z ** 2.0 - other.OriginalCenter.Z ** 2.0 -
                        this.OriginalRadius ** 2.0 + other. OriginalRadius ** 2.0

                // Intersection center
                let t = (this.OriginalCenter.X * A + this.OriginalCenter.Y * B + this.OriginalCenter.Z * C + D) /
                        (A * (this.OriginalCenter.X - other.OriginalCenter.X) +
                         B * (this.OriginalCenter.Y - other.OriginalCenter.Y) +
                         C * (this.OriginalCenter.Z - other.OriginalCenter.Z))
                let x = this.OriginalCenter.X + t * (other.OriginalCenter.X - this.OriginalCenter.X)
                let y = this.OriginalCenter.Y + t * (other.OriginalCenter.Y - this.OriginalCenter.Y)
                let z = this.OriginalCenter.Z + t * (other.OriginalCenter.Z - this.OriginalCenter.Z)
                let intersectionCenter: Point = { X = x; Y = y; Z = z }

                // Intersection
                let alpha = Math.Acos((this.OriginalRadius ** 2.0 + dist ** 2.0 - other.OriginalRadius ** 2.0) /
                                      (2.0 * this.OriginalRadius * dist))
                let R = this.OriginalRadius * Math.Sin alpha

                match R with
                | 0.0 -> IntersectionPoint intersectionCenter
                | _ ->
                    let v = this.OriginalCenter.FindVector other.OriginalCenter
                    IntersectionCircle (intersectionCenter, R, v)

    type ViewBox = float * float * float * float

    type Depiction =
        | Filled
        | BallAndStick

    let origin: Point = { X = 0.0; Y = 0.0; Z = 0.0 }

    let filterAtoms (atomType: Atom) (atoms: AtomInfo array) : AtomInfo array =
        Array.filter(fun (atom: AtomInfo) -> atom.Type <> atomType) atoms

    let physicalProjection cameraPerpendicular cameraHorizon cameraForward (pov: Point) (p: Point)  : Point =
        let pointVector = pov.FindVector(p)
        { X = pointVector.ProjectVector cameraPerpendicular
          Y = pointVector.ProjectVector cameraHorizon
          Z = pointVector.ProjectVector cameraForward }

    let perspectiveProjection (focalLength: float) (p: Point) : Point =
        let scaleFactor = focalLength / p.Z
        { X = p.X * scaleFactor
          Y = p.Y * scaleFactor
          Z = p.Z * scaleFactor }

    let project cameraPerpendicular cameraHorizon cameraForward (pov: Point) focalLength (p: Point) : Point =
        p |> physicalProjection cameraPerpendicular cameraHorizon cameraForward pov |> perspectiveProjection focalLength

    let createAtom (index: int) (atomType: Atom) (center: Point) (radius: Radius) : AtomInfo =
        { Index = index
          Type = atomType
          OriginalCenter = center
          OriginalRadius = radius
          ProjectedCenter = center
          ProjectedRadius = radius }