namespace CineMol.Types


module Fundamentals =

    /// <summary>
    /// Index indicates the indexation number of an object.
    /// </summary>
    type Index = Index of int

    /// <summary>
    /// Radius defines the radius of a circle or sphere.
    /// </summary>
    type Radius = Radius of float


module Style =

    /// <summary>
    /// Color describes a color in RGB int values or Hex string.
    /// </summary>
    type Color = Color of int * int * int
        with
        member this.ToHex =
            let (Color (r, g, b)) = this
            sprintf "#%02x%02x%02x" <| r, g, b


module Geometry =

    open System

    open Fundamentals

    /// <summary>
    /// Vector2D resembles a vector in two-dimensional Euclidean space.
    /// </summary>
    type Vector2D = { X: float; Y: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X;  Y = p1.Y + p2.Y }
        static member (-) (p1, p2) = { X = p1.X - p2.X; Y = p1.Y - p2.Y }
        member this.Mul v = { X = this.X * v; Y = this.Y * v }
        member this.Div v = { X = this.X / v; Y = this.Y / v }
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v }
        member this.Dot other = this.X * other.X + this.Y + other.Y
        member this.Mag = this.SumOfSquares |> Math.Sqrt
        member this.Sum = this.X + this.Y
        member this.SumOfSquares = (this.Pow 2.0).Sum
        member this.Norm = this.Mul (if this.Mag = 0.0 then infinity else 1.0 / this.Mag)

    /// <summary>
    /// Vector3D resembles a vector in three-dimensional Euclidean space.
    /// </summary>
    type Vector3D = { X: float; Y: float; Z: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y; Z = p1.Z + p2.Z }
        static member (-) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y; Z = p1.Z + p2.Z }
        member this.Mul v = { X = this.X * v; Y = this.Y * v; Z = this.Z * v }
        member this.Div v = { X = this.X / v; Y = this.Y / v; Z = this.Z / v }
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v; Z = this.Z ** v }
        member this.Dot other = this.X * other.X + this.Y * other.Y + this.Z * other.Z
        member this.Mag = this.SumOfSquares |> Math.Sqrt
        member this.Sum = this.X + this.Y + this.Z
        member this.SumOfSquares = (this.Pow 2.0).Sum
        member this.Norm = this.Mul (if this.Mag = 0.0 then infinity else 1.0 / this.Mag)
        member this.Cross other =
            { X = this.Y * other.Z - this.Z * other.Y
              Y = this.Z * other.X - this.X * other.Z
              Z = this.X * other.Y - this.Y * other.X }

    /// <summary>
    /// Point2D resembles a point in two-dimensional Euclidean space.
    /// </summary>
    type Point2D = { X: float; Y: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y }
        static member (-) (p1, p2) = { X = p1.X - p2.X; Y = p1.Y - p2.Y }
        member this.Mul v = { X = this.X * v; Y = this.X * v }
        member this.Div v = { X = this.X / v; Y = this.Y / v }
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v }
        member this.Sum = this.X + this.Y
        member this.Distance other = ((this - other).Pow 2.0).Sum |> Math.Sqrt
        member this.Midpoint other = (this + other).Div 2.0
        member this.FindVector other = other - this

    /// <summary>
    /// Point3D resembles a point in three-dimensional Euclidean space.
    /// </summary>
    type Point3D = { X: float; Y: float; Z: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y; Z = p1.Z + p2.Z }
        static member (-) (p1, p2) = { X = p1.X - p2.X; Y = p1.Y - p2.Y; Z = p1.Z - p2.Z }
        member this.Mul v = { X = this.X ** v; Y = this.Y ** v; Z = this.Z ** v }
        member this.Div v = { X = this.X / v; Y = this.Y / v; Z = this.Z / v }
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v; Z = this.Z ** v }
        member this.Sum = this.X + this.Y + this.Z
        member this.Distance other = ((this - other).Pow 2.0).Sum |> Math.Sqrt
        member this.Midpoint other = (this + other).Div 2.0
        member this.FindVector other = other - this
        member p.Rotate (axis: Axis) rad = axis.RotationMatrix(p, rad)
        member this.ToPoint2D () = { X = this.X; Y = this.Y }

    /// <summary>
    /// Axis describes a plane in a three-dimensional Euclidean space.
    /// </summary>
    and Axis = | X | Y | Z
        with
        member this.RotationMatrix =
            match this with
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

    /// <summary>
    /// Definition for a line.
    /// </summary>
    type Line = Line of Point2D * Point2D

    /// <summary>
    /// Definition for a circle.
    /// </summary>
    type Circle = Circle of Point2D * Radius

    /// <summary>
    /// Definition for a quadrangle.
    /// </summary>
    type Quadrangle = Quadrangle of Point2D * Point2D * Point2D * Point2D

    /// <summary>
    /// Definition for a sphere.
    /// </summary>
    type Sphere = Sphere of Point3D * Radius

    /// <summary>
    /// Definition for a cylinder.
    /// </summary>
    type Cylinder = Cylinder of Line * Radius


module Chem =

    open Fundamentals
    open Style
    open Geometry

    /// <summary>
    /// AtomType describes the atomic number of an atom.
    /// </summary>
    type AtomType =
        | H                                                                                  | He
        | Li | Be                                                   | B  | C  | N  | O  | F  | Ne
        | Na | Mg                                                   | Al | Si | P  | S  | Cl | Ar
        | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn | Ga | Ge | As | Se | Br | Kr
        | Rb | Sr | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | Sb | Te | I  | Xe
        | Cs | Ba | Lu | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg | Tl | Pb | Bi | Po | At | Rn
        | Fr | Ra | Lr | Rf | Db | Sg | Bh | Hs | Mt | Ds | Rg | Cn | Nh | Fl | Mc | Lv | Ts | Og

    /// <summary>
    /// AtomInfo records all information on the atom identity and styling.
    /// </summary>
    type AtomInfo =
        { Index: Index
          Type: AtomType
          Color: Color option }

    /// <summary>
    /// BondInfo records all information on the bond identity and styling.
    /// </summary>
    type BondInfo =
        { Index: Index
          BeginAtom: Index
          EndAtom: Index
          Color: Color option }

    /// <summary>
    /// Atom describes an atom in two-dimensional or three-dimensional
    /// Euclidean space.
    /// </summary>
    type Atom =
        | Atom2D of AtomInfo * Circle
        | Atom3D of AtomInfo * Sphere

    /// <summary>
    /// Bond describes a bond in two-dimensional or three-dimensional
    /// Euclidean space.
    /// </summary>
    type Bond =
        | Bond2D of BondInfo * Line
        | Bond3D of BondInfo * Cylinder


module Svg =

    open Geometry

    /// <summary>
    /// ViewBox defines the boundaries of the SVG viewbox.
    /// </summary>
    type ViewBox = { MinX: float; MinY: float; Width: float; Height: float }

    /// <summary>
    /// Definition is a collection of supported SVG definition embeddings.
    /// SVG definitions are elements that can be reused in inside an SVG image.
    /// </summary>
    type Definition = | RadialGradient | LinearGradient

    /// <summary>
    /// Shape is a collection of supported shapes to draw in
    /// two-dimensional Euclidean space as SVG XML objects.
    /// </summary>
    type Shape =
        | Line of Line
        | Circle of Circle
        | Quadrangle of Quadrangle

    /// <summary>
    /// Header describes the SVG ID and the SVG viewbox.
    /// </summary>
    type Header = Header of ID: string * ViewBox

    /// <summary>
    /// SVG encapsulates all individual elements in the SVG image.
    /// </summary>
    type SVG =
        { Header: Header
          Definitions: Definition list
          Objects: Shape list }
