namespace CineMol.Types

open System
open CineMol.Helpers
    
/// <summary>
/// Color describes a color in RGB int values or Hex string.
/// </summary>
type Color = Color of int * int * int
    with
    override this.ToString () =
        let (Color (r, g, b)) = this 
        $"rgb(%i{r},%i{g},%i{b})"
    
    member this.ToHex =
        let (Color (r, g, b)) = this
        $"#{r:x2}{g:x2}{b:x2}"

    member this.Diffuse alpha =
        let alpha = clamp 0.0 1.0 alpha
        let (Color (r, g, b)) = this
        let diffuseChannel c = (float c) * alpha |> int
        ( diffuseChannel r,
          diffuseChannel g,
          diffuseChannel b ) |> Color
        
/// <summary>
/// Model styles.
/// </summary>
type ModelStyle = | SpaceFilling 

module Geometry =

    /// <summary>
    /// Point2D resembles a point in two-dimensional Euclidean space.
    /// </summary>
    type Point2D = { X: float; Y: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y }
        
        static member (-) (p1, p2) = { X = p1.X - p2.X; Y = p1.Y - p2.Y }
        
        static member (*) (p1, p2) = { X = p1.X * p2.X; Y = p1.Y * p2.Y }
        
        member this.Add v = { X = this.X + v; Y = this.Y + v }
        
        member this.Mul v = { X = this.X * v; Y = this.X * v }
        
        member this.Div v = { X = this.X / v; Y = this.Y / v }
        
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v }
        
        member this.Sum () = this.X + this.Y
        
        member this.Dist other = ((this - other).Pow 2.0).Sum() |> Math.Sqrt
        
        member this.Midpoint other = (this + other).Div 2.0
        
        member this.CreateVector (other: Point2D) : Vector2D = { X = other.X - this.X; Y = other.Y - this.Y }

        static member Centroid (ps: Point2D list) =
            ps |> List.fold (fun pSum p -> pSum + p) { X = 0.0; Y = 0.0 } |> (fun p -> p.Div (float ps.Length))
            
    and Vector2D = { X: float; Y: float }
        with
        member this.Cross (other: Vector2D) = this.X * other.Y - this.Y * other.X
        
    /// <summary>
    /// Point3D resembles a point in three-dimensional Euclidean space.
    /// </summary>
    and Point3D = { X: float; Y: float; Z: float }
        with
        static member (+) (p1, p2) = { X = p1.X + p2.X; Y = p1.Y + p2.Y; Z = p1.Z + p2.Z }
        
        static member (-) (p1, p2) = { X = p1.X - p2.X; Y = p1.Y - p2.Y; Z = p1.Z - p2.Z }
        
        static member (*) (p1, p2) = { X = p1.X * p2.X; Y = p1.Y * p2.Y; Z = p1.Z * p2.Z }
        
        member this.Add v = { X = this.X + v; Y = this.Y + v; Z = this.Z + v }
        
        member this.Mul v = { X = this.X * v; Y = this.Y * v; Z = this.Z * v }
        
        member this.Div v = { X = this.X / v; Y = this.Y / v; Z = this.Z / v }
        
        member this.Pow v = { X = this.X ** v; Y = this.Y ** v; Z = this.Z ** v }
        
        member this.Sum () = this.X + this.Y + this.Z
        
        member this.Dist other = ((this - other).Pow 2.0).Sum() |> Math.Sqrt
        
        member this.Midpoint other = (this + other).Div 2.0

        member p.Rotate axis rad =
            match axis with
            | X ->
               { X = p.X
                 Y = p.Y * Math.Cos(rad) - p.Z * Math.Sin(rad)
                 Z = p.Y * Math.Sin(rad) + p.Z * Math.Cos(rad) }
            | Y ->
               { X = p.X * Math.Cos(rad) + p.Z * Math.Sin(rad)
                 Y = p.Y
                 Z = p.Z * Math.Cos(rad) - p.X * Math.Sin(rad) }
            | Z ->
               { X = p.X * Math.Cos(rad) - p.Y * Math.Sin(rad)
                 Y = p.X * Math.Sin(rad) + p.Y * Math.Cos(rad)
                 Z = p.Z }
               
        member this.ToPoint2D () = { X = this.X; Y = this.Y }
        
        static member Centroid (ps: Point3D list) =
            ps |> List.fold (fun pSum p -> pSum + p) { X = 0.0; Y = 0.0; Z = 0.0 } |> (fun p -> p.Div (float ps.Length))
    
    /// <summary>
    /// Axis describes a plane in a three-dimensional Euclidean space.
    /// </summary>
    and Axis = | X | Y | Z
        with
        static member Origin () = { X = 0.0; Y = 0.0; Z = 0.0 }
    
    /// <summary>
    /// Definition for a circle in two-dimensional Euclidean space.
    /// </summary>
    and Circle =
        { Center : Point2D
          Radius : float }
    
    /// <summary>
    /// Definition for a sphere.
    /// </summary>
    and Sphere =
        { Center: Point3D
          Radius : float }
        with
        member this.Encloses (o: Sphere) : bool =
            let dist = this.Center.Dist o.Center
            dist + o.Radius <= this.Radius
            
        member this.Intersects (o: Sphere) : bool =
            let dist = this.Center.Dist o.Center
            dist <= this.Radius + o.Radius
            
        member this.PointsOnSphere (resolution : int) : Point3D list =
            let N = float resolution 
            
            // Calculate polar angles.
            let numPointsPhi = N / 2.0 |> int
            let phis =
                [ for i in [ 1 .. 1 .. numPointsPhi ] do
                    yield ((float i) / (( float numPointsPhi) - 1.0)) * Math.PI ]
            
            // Calculate azimuthal angles.
            let numPointsTheta = N / 2.0 |> int
            let thetas =
                [ for i in [ 0 .. 1 .. numPointsTheta ] do
                    yield ((float i) / (float numPointsTheta)) * 2.0 * Math.PI ]
            
            // Calculate points on sphere.
            let points = 
                [
                    for phi in phis do
                        
                        // Only calculate points on half of sphere we can see over the z-axis.
                        let z = this.Center.Z + this.Radius * Math.Cos(phi)
                        match z with
                        | z when z < this.Center.Z -> () // We can never see this part of the sphere.
                        | _ ->
                            // let points : Point3D list =
                            //     [
                            //         for theta in thetas do
                            //             let x = this.Center.X + this.Radius * Math.Sin(phi) * Math.Cos(theta)
                            //             let y = this.Center.Y + this.Radius * Math.Sin(phi) * Math.Sin(theta)
                            //             yield { X = x; Y = y; Z = z } 
                            //     ]
                            //     |> List.rev // Point furthest away from us first.
                            //
                            // match points with
                            // | [] -> ()
                            // | _ -> yield points
                            
                            for theta in thetas do 
                                let x = this.Center.X + this.Radius * Math.Sin(phi) * Math.Cos(theta)
                                let y = this.Center.Y + this.Radius * Math.Sin(phi) * Math.Sin(theta)
                                yield { X = x; Y = y; Z = z }
                ]
                
            // Transpose the list of lists.
            // quadrants |> List.transpose
            
            points 

module Chem =

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
        | Fr | Ra
        with
        static member FromString (atomString: string) =
            match atomString with
            | "H"  -> Some H  | "He" -> Some He | "Li" -> Some Li | "Be" -> Some Be | "B"  -> Some B 
            | "C"  -> Some C  | "N"  -> Some N  | "O"  -> Some O  | "F"  -> Some F  | "Ne" -> Some Ne
            | "Na" -> Some Na | "Mg" -> Some Mg | "Al" -> Some Al | "Si" -> Some Si | "P"  -> Some P
            | "S"  -> Some S  | "Cl" -> Some Cl | "Ar" -> Some Ar | "K"  -> Some K  | "Ca" -> Some Ca
            | "Sc" -> Some Sc | "Ti" -> Some Ti | "V"  -> Some V  | "Cr" -> Some Cr | "Mn" -> Some Mn
            | "Fe" -> Some Fe | "Co" -> Some Co | "Ni" -> Some Ni | "Cu" -> Some Cu | "Zn" -> Some Zn
            | "Ga" -> Some Ga | "Ge" -> Some Ge | "As" -> Some As | "Se" -> Some Se | "Br" -> Some Br
            | "Kr" -> Some Kr | "Rb" -> Some Rb | "Sr" -> Some Sr | "Y"  -> Some Y  | "Zr" -> Some Zr
            | "Nb" -> Some Nb | "Mo" -> Some Mo | "Tc" -> Some Tc | "Ru" -> Some Ru | "Rh" -> Some Rh
            | "Pd" -> Some Pd | "Ag" -> Some Ag | "Cd" -> Some Cd | "In" -> Some In | "Sn" -> Some Sn
            | "Sb" -> Some Sb | "Te" -> Some Te | "I"  -> Some I  | "Xe" -> Some Xe | "Cs" -> Some Cs
            | "Ba" -> Some Ba | "Lu" -> Some Lu | "Hf" -> Some Hf | "Ta" -> Some Ta | "W"  -> Some W
            | "Re" -> Some Re | "Os" -> Some Os | "Ir" -> Some Ir | "Pt" -> Some Pt | "Au" -> Some Au
            | "Hg" -> Some Hg | "Tl" -> Some Tl | "Pb" -> Some Pb | "Bi" -> Some Bi | "Po" -> Some Po
            | "At" -> Some At | "Rn" -> Some Rn | "Fr" -> Some Fr | "Ra" -> Some Ra | _    -> None 
    
    /// <summary>
    /// BondType describes the type of bond between two atoms.
    /// </summary>
    type BondType = | Single | Double | Triple | Aromatic
        with 
        static member FromString (bondString: string) =
            match bondString with
            | "1" | "SINGLE"   | "Single"   -> Some Single
            | "2" | "DOUBLE"   | "Double"   -> Some Double
            | "3" | "TRIPLE"   | "Triple"   -> Some Triple
            | "4" | "AROMATIC" | "Aromatic" -> Some Aromatic
            | _                             -> None 
    
    /// <summary>
    /// BondInfo records all information on the bond identity and styling.
    /// </summary>
    type BondInfo =
        { BeginAtomIndex: int
          EndAtomIndex: int
          Type: BondType 
          Color: Color option }
    
    /// <summary>
    /// Atom3D describes an atom in three-dimensional Euclidean space.
    /// </summary>
    type Atom =
        { Index : int
          Type : AtomType
          Color : Color
          Position : Point3D
          Radius : float }
    
    /// <summary>
    /// Bond describes a bond between two Atoms in two-dimensional or three-dimensional Euclidean space.
    /// </summary>
    type Bond = Bond of BondInfo
        with
        member this.Unwrap () =
            let (Bond info) = this
            info 
    
    /// <summary>
    /// Molecule describes a molecule, which contains of Atoms and Bonds.
    /// </summary>
    type Molecule = { Atoms: Atom list; Bonds: Bond list }
        with
        member this.AdjustForCentroid () =
            let centroid = this.Atoms |> List.map (fun atom -> atom.Position) |> Point3D.Centroid
            let adjustedAtoms = this.Atoms |> List.map (fun atom -> { atom with Position = atom.Position - centroid })
            { this with Atoms = adjustedAtoms }  

module Svg =
    
    open Geometry
    
    /// <summary>
    /// ViewBox defines the boundaries of the SVG view box.
    /// </summary>
    type ViewBox = { MinX: float; MinY: float; Width: float; Height: float }
        with
        override this.ToString () =
            $"viewBox=\"%.3f{this.MinX} %.3f{this.MinY} %.3f{this.Width} %.3f{this.Height}\""

    /// <summary>
    /// Shape is a collection of supported shapes to draw in two-dimensional Euclidean space as SVG XML objects.
    /// </summary>
    type Shape =
        | Circle of Index * Color * Circle
        | Polygon of Index * Color * Point2D list
        with
        member this.ToSvg () =
            match this with
            | Circle (index, color, circle) ->
                let x = circle.Center.X
                let y = circle.Center.Y
                $"<circle class=\"item-{index}\" cx=\"%.3f{x}\" cy=\"%.3f{y}\" r=\"%.3f{circle.Radius}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.1\"/>"

            | Polygon (index, color, points) ->
                let pointsStr = points |> List.map (fun p -> $"%.3f{p.X},%.3f{p.Y}") |> String.concat " "
                $"<polygon class=\"item-{index}\" points=\"{pointsStr}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.1\"/>"
        
    and Index = int 
                
    /// <summary>
    /// Header describes the SVG ID and the SVG view box.
    /// </summary>
    type Header = Header of version: float * encoding: string 
        with
        static member New () = Header (1.0, "UTF-8")
        override this.ToString () : string =
            let (Header (version, encoding)) = this
            $"<?xml version=\"%.1f{version}\" encoding=\"{encoding}\"?>"
    
    /// <summary>
    /// SVG encapsulates all individual elements in the SVG image.
    /// </summary>
    type SVG = { Header: Header; ID: string; ViewBox: ViewBox; Objects: Shape list }
        with
        override this.ToString () =            
            // Concatenate definitions, objects, and header strings. 
            this.Header.ToString() + this.Body() 
            
        member this.Body() =            
            // Convert all objects to a single string.
            let objs = this.Objects |> List.map (fun x -> x.ToSvg()) |> String.concat "\n"
            
            // Combine items into SVG body.
            $"<svg id=\"{this.ID}\" xmlns=\"http://www.w3.org/2000/svg\" {this.ViewBox.ToString()}>\n{objs}\n</svg>"
            
module Drawing =
    
    open Svg
    
    /// <summary>
    /// Drawing options.
    /// </summary>
    type DrawingOptions =
        { ViewBox: ViewBox option
          Style: ModelStyle
          Resolution: int }
        with
        static member New () =
            { ViewBox = None
              Style = SpaceFilling
              Resolution = 50 }