module Client.CineMol.Types

open System
open Styles

type Zoom = { Ratio: float }
    with
    static member init = { Ratio = 1.0 }

type Rotation = { AxisX: float; AxisY: float; AxisZ: float }
    with
    static member init =
        { AxisX = 0.0; AxisY = 0.0; AxisZ = 0.0 }

type Axis = | X | Y | Z
    with
    member x.RotationMatrix : Point3D * float -> Point3D =
        match x with
        | X ->
            (fun (p: Point3D, rad: float) ->
                { X = p.X
                  Y = p.Y * Math.Cos(rad) - p.Z * Math.Sin(rad)
                  Z = p.Y * Math.Sin(rad) + p.Z * Math.Cos(rad) })
        | Y ->
            (fun (p: Point3D, rad: float) ->
                { X = p.X * Math.Cos(rad) + p.Z * Math.Sin(rad)
                  Y = p.Y
                  Z = p.Z * Math.Cos(rad) - p.X * Math.Sin(rad) })
        | Z ->
            (fun (p: Point3D, rad: float) ->
                { X = p.X * Math.Cos(rad) - p.Y * Math.Sin(rad)
                  Y = p.X * Math.Sin(rad) + p.Y * Math.Cos(rad)
                  Z = p.Z })

and Point2D = { X: float; Y: float }
    with
    static member (-) (p1: Point2D, p2: Point2D) : Point2D =
        { X = p1.X - p2.X; Y = p1.Y - p2.Y }

    static member Pow (p: Point2D) (d: float) : Point2D =
        { X = p.X ** d; Y = p.Y ** d }

    static member Sum (p: Point2D) : float =
        p.X + p.Y

    member p1.Distance (p2: Point2D) : float =
        Math.Sqrt(Point2D.Sum(Point2D.Pow (p1 - p2) 2.0))

    member p1.Centroid (p2: Point2D) : Point2D =
        { X = (p1.X + p2.X) / 2.0
          Y = (p1.Y + p2.X) / 2.0 }

    member p1.FindVector (p2: Point2D) : Vector2D =
        { X = p2.X - p1.X; Y = p2.Y - p1.Y }

and Vector2D = { X: float; Y: float }
    with
    member u.SumOfSquares : float =
        u.X ** 2.0 + u.Y ** 2.0

    member u.Magnitude : float =
        Math.Sqrt(u.SumOfSquares)

    static member (*) (k, v: Vector2D) =
        { X = k * v.X; Y = k * v.Y }

    static member (+) (v1: Vector2D, v2: Vector2D) =
        { X = v1.X + v2.X; Y = v1.Y + v2.Y }

    member u.Norm : Vector2D =
        let mag = u.Magnitude
        let div = if mag = 0.0 then infinity else 1.0 / mag
        div * u

    member u.Dot (v: Vector2D) : float = u.X * v.X + u.Y * v.Y

and Point3D = { X: float; Y: float; Z: float }
    with
    static member (-) (p1: Point3D, p2: Point3D) : Point3D =
        { X = p1.X - p2.X; Y = p1.Y - p2.Y; Z = p1.Z - p2.Z }

    static member Pow (p: Point3D) (d: float) : Point3D =
        { X = p.X ** d; Y = p.Y ** d; Z = p.Z ** d }

    static member Sum (p: Point3D) : float =
        p.X + p.Y + p.Z

    member p1.Distance (p2: Point3D) : float =
        Math.Sqrt(Point3D.Sum(Point3D.Pow (p1 - p2) 2.0))

    member p1.Centroid (p2: Point3D) : Point3D =
        { X = (p1.X + p2.X) / 2.0
          Y = (p1.Y + p2.X) / 2.0
          Z = (p1.Z + p2.Z) / 2.0 }

    member p.Rotate (axis: Axis) (rad: float) : Point3D =
        axis.RotationMatrix(p, rad)

    member p1.FindVector (p2: Point3D) : Vector3D =
        { X = p2.X - p1.X; Y = p2.Y - p1.Y; Z = p2.Z - p1.Z }

and Vector3D = { X: float; Y: float; Z: float }
    with
    member u.SumOfSquares : float =
        u.X ** 2.0 + u.Y ** 2.0 + u.Z ** 2.0

    member u.Magnitude : float =
        Math.Sqrt(u.SumOfSquares)

    static member (*) (k, v: Vector3D) =
        { X = k * v.X; Y = k * v.Y; Z = k * v.Z }

    static member (+) (v1: Vector3D, v2: Vector3D) =
        { X = v1.X + v2.X; Y = v1.Y + v2.Y; Z = v1.Z + v2.Z }

    member u.Norm : Vector3D =
        let mag = u.Magnitude
        let div = if mag = 0.0 then infinity else 1.0 / mag
        div * u

    member u.Dot (v: Vector3D) : float = u.X * v.X + u.Y * v.Y + u.Z * v.Z

    member u.Cross (v: Vector3D) : Vector3D =
        { X = u.Y * v.Z - u.Z * v.Y
          Y = u.Z * v.X - u.X * v.Z
          Z = u.X * v.Y - u.Y * v.X }

    member x.ProjectVector (v: Vector3D) : float = (v.Dot x) / v.Magnitude

type Index = int

type SphereSphereIntersection =
    | Eclipsed
    | NoIntersection
    | IntersectionPoint of Point3D
    | IntersectionCircle of Point3D * Radius * Vector3D

type Clipping = { Line: Line }
and Line = Point2D * Point2D

type AtomInfo =
    { Index: Index
      AtomType: AtomType
      Center: Point3D
      Radius: Radius }
    with
    member x.Rotate (axis: Axis) (rad: float) : AtomInfo =
        { x with Center = x.Center.Rotate axis rad }

    member this.Intersect (other: AtomInfo) : SphereSphereIntersection =
        let dist = this.Center.Distance other.Center

        match dist with
        | d when d >= (this.Radius + other.Radius) || (d = 0.0 && this.Radius = other.Radius)
            -> NoIntersection

        | d when (d + this.Radius) < other.Radius
            -> Eclipsed

        | _ ->
            // Intersection plane
            let A = 2.0 * (other.Center.X - this.Center.X)
            let B = 2.0 * (other.Center.Y - this.Center.Y)
            let C = 2.0 * (other.Center.Z - this.Center.Z)
            let D = this.Center.X ** 2.0 - other.Center.X ** 2.0 + this.Center.Y ** 2.0 - other.Center.Y ** 2.0 +
                    this.Center.Z ** 2.0 - other.Center.Z ** 2.0 - this.Radius ** 2.0 + other.Radius ** 2.0

            // Intersection center
            let t = (this.Center.X * A + this.Center.Y * B + this.Center.Z * C + D) /
                    (A * (this.Center.X - other.Center.X) + B * (this.Center.Y - other.Center.Y) + C * (this.Center.Z - other.Center.Z))
            let x = this.Center.X + t * (other.Center.X - this.Center.X)
            let y = this.Center.Y + t * (other.Center.Y - this.Center.Y)
            let z = this.Center.Z + t * (other.Center.Z - this.Center.Z)
            let intersectionCenter: Point3D = { X = x; Y = y; Z = z }

            // Intersection
            let x = (this.Radius ** 2.0 + dist ** 2.0 - other.Radius ** 2.0) / (2.0 * this.Radius * dist)
            if x < 1.0 then
                let alpha = Math.Acos(x)
                let R = this.Radius * Math.Sin alpha
                match R with
                | 0.0 -> IntersectionPoint intersectionCenter
                | _ ->
                    let v = this.Center.FindVector other.Center
                    IntersectionCircle (intersectionCenter, R, v)
            else
                NoIntersection

type ProjectedAtomInfo =
    { Index: int
      AtomType: AtomType
      Center: Point2D
      Radius: Radius
      Clipping: Clipping list }


type BondInfo =
    { Index: Index
      Start: Index
      End: Index
      BondType: BondType
      Scaling: float }
and BondType = | Single | Double | Triple | Aromatic | Unknown

let createAtom (index: int) (atomType: AtomType) (c: Point3D) (r: Radius) : AtomInfo =
    { Index = index; AtomType = atomType; Center = c; Radius = r }

type Molecule = { Atoms: AtomInfo[]; Bonds: BondInfo[] }
type ProjectedMolecule = { Atoms: ProjectedAtomInfo[]; Bonds: BondInfo[] }

type ViewBox = float * float * float * float

type Depiction = | Filled | BallAndStick | Wire

let origin: Point3D = { X = 0.0; Y = 0.0; Z = 0.0 }

let physicalProjection
    cameraPerpendicular
    cameraHorizon
    cameraForward
    (pov: Point3D)
    (p: Point3D) : Point3D =
    let pointVector: Vector3D = pov.FindVector p
    { X = pointVector.ProjectVector cameraPerpendicular
      Y = pointVector.ProjectVector cameraHorizon
      Z = pointVector.ProjectVector cameraForward }

let perspectiveProjection (focalLength: float) (p: Point3D) : Point3D =
    let scaleFactor = focalLength / p.Z
    { X = p.X * scaleFactor
      Y = p.Y * scaleFactor
      Z = p.Z * scaleFactor }

let project
    cameraPerpendicular
    cameraHorizon
    cameraForward
    (pov: Point3D)
    focalLength
    (p: Point3D)
    : Point3D =
    p
    |> physicalProjection cameraPerpendicular cameraHorizon cameraForward pov
    |> perspectiveProjection focalLength
