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
type ModelStyle =
    | SpaceFilling
    | BallAndStick
    | WireFrame
    
/// <summary>
/// Art styles.
/// </summary>
type ArtStyle =
    | Cartoon
    | Glossy 

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
        
        member this.VectorTo other : Vector3D = { X = other.X - this.X; Y = other.Y - this.Y; Z = other.Z - this.Z }
        
        member this.MoveTowards other (distance: float) =
            let dist = this.Dist other
            let vector = this.VectorTo other
            match dist - distance with
            | newDist when newDist <= 0.0 -> other
            | newDist ->
                let norm = vector.Normalize ()
                this + { X = norm.X * newDist; Y = norm.Y * newDist; Z = norm.Z * newDist }
            
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
    
    and Vector3D = { X: float; Y: float; Z: float }
        with
        member this.Normalize () =
            let length = Math.Sqrt (this.X * this.X + this.Y * this.Y + this.Z * this.Z)
            { X = this.X / length; Y = this.Y / length; Z = this.Z / length }
    
    /// <summary>
    /// Axis describes a plane in a three-dimensional Euclidean space.
    /// </summary>
    and Axis = | X | Y | Z
        with
        static member Origin () : Point3D = { X = 0.0; Y = 0.0; Z = 0.0 }
    
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
                                let point : Point3D = { X = x; Y = y; Z = z }
                                yield point 
                ]
                
            // Transpose the list of lists.
            // quadrants |> List.transpose
            
            points
            
    /// <summary>
    /// Definition for a cylinder.
    /// </summary>
    type Cylinder =
        { Start : Point3D
          End : Point3D
          Radius : float }
        with
        member this.IsInside (point: Point3D) : bool =
            let cylinderDirection = this.End - this.Start
            let pointToStart = point - this.Start
            let projection =
                let nominator = pointToStart |> fun p -> p.X * cylinderDirection.X + p.Y * cylinderDirection.Y + p.Z * cylinderDirection.Z
                let denominator = cylinderDirection |> fun c -> c.X * c.X + c.Y * c.Y + c.Z * c.Z
                nominator / denominator
                
            if projection < 0.0 then false // The point is behind the cylinder's start point
            elif projection > 1.0 then false // The point is beyond the cylinder's end point
            else
                let closestPointOnAxis = this.Start.Add(projection) * cylinderDirection
                let distanceSquared = (point - closestPointOnAxis) |> fun d -> d.X * d.X + d.Y * d.Y + d.Z * d.Z
                let radiusSquared = this.Radius * this.Radius
                distanceSquared <= radiusSquared
                
        member this.PointsOnCylinder (resolution : int) : Point3D list =
            let N = resolution
            
            let latitudeDivisions = N / 2
            let longitudeDivisions = N
            
            [
                for lat = 0 to latitudeDivisions do
                    let theta = float lat * Math.PI / float latitudeDivisions
                    let sinTheta = sin theta
                    let cosTheta = cos theta

                    for lon = 0 to longitudeDivisions do
                        let phi = float lon * 2.0 * Math.PI / float longitudeDivisions
                        let sinPhi = sin phi
                        let cosPhi = cos phi

                        let x = this.Radius * sinTheta * cosPhi
                        let y = this.Radius * sinTheta * sinPhi
                        let z = this.Radius * cosTheta

                        let point : Point3D = { X = x; Y = y; Z = z }
                        yield point 
            ]
            
    /// <summary>
    /// Calculate the centroid of a list of points.
    /// </summary>
    let calcCentroid (points : Point2D list) : Point2D =
        let start : Point2D = { X = 0.0; Y = 0.0 }
        let sum = points |> List.fold (fun pSum p -> pSum + p) start
        let count = float points.Length
        sum.Div count 

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
    /// Atom3D describes an atom in three-dimensional Euclidean space.
    /// </summary>
    type Atom =
        { Index : int
          Type : AtomType
          Color : Color
          Opacity : float 
          Position : Point3D
          Radius : float }
    
    /// <summary>
    /// Bond describes a bond between two Atoms in two-dimensional or three-dimensional Euclidean space.
    /// </summary>
    type Bond =
        { BeginIndex : int
          EndIndex : int 
          Type : BondType
          BeginAtomIndex : int
          EndAtomIndex : int
          Opacity : float option 
          Color : Color option
          Radius : float }
    
    /// <summary>
    /// Molecule describes a molecule, which contains of Atoms and Bonds.
    /// </summary>
    type Molecule = { Atoms: Atom list; Bonds: Bond list }
        with
        member this.AdjustForCentroid () =
            let centroid = this.Atoms |> List.map (fun atom -> atom.Position) |> Point3D.Centroid
            let adjustedAtoms = this.Atoms |> List.map (fun atom -> { atom with Position = atom.Position - centroid })
            { this with Atoms = adjustedAtoms }
            
        member this.GetAtom (atomIndex: int) : Atom option =
            try this.Atoms |> List.find (fun atom -> atom.Index = atomIndex) |> Some 
            with _ -> None
            
        member this.GetBonds (includeHydrogenAtoms : bool) (atomIndex: int) : Bond list =
            this.Bonds
            |> List.filter (fun bond -> bond.BeginAtomIndex = atomIndex || bond.EndAtomIndex = atomIndex)
            |> List.filter (fun bond ->
                let beginAtom = this.GetAtom bond.BeginAtomIndex
                let endAtom = this.GetAtom bond.EndAtomIndex
                match beginAtom, endAtom with
                | Some beginAtom, Some endAtom ->
                    if includeHydrogenAtoms then true
                    else beginAtom.Type <> H && endAtom.Type <> H
                | _ -> false)

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
    /// SVG fills.
    /// </summary>
    type Fill =
        | RadialGradient of int * Point2D * float * Color
        | LinearGradient of int * Point2D * Point2D * Color
        with
        override this.ToString () =
            match this with
            | LinearGradient (index, start, stop, color) ->
                let x1 = start.X
                let y1 = start.Y
                let x2 = stop.X
                let y2 = stop.Y
                let startColor = color
                let stopColor = color.Diffuse 0.5
                $"<linearGradient id=\"item-{index}\" x1=\"%.3f{x1}\" x2=\"%.3f{x2}\" y1=\"%.3f{y1}\" y2=\"%.3f{y2}\" gradientUnits=\"userSpaceOnUse\" spreadMethod=\"reflect\"><stop offset=\"0.00\" stop-color=\"{startColor}\"/><stop offset=\"1.00\" stop-color=\"{stopColor}\"/></linearGradient>"
            
            | RadialGradient (index, center, radius, color) ->
                let cx = center.X
                let cy = center.Y
                let r = radius * 1.5 // Make the gradient a bit larger than the object.
                let startColor = color 
                let stopColor = color.Diffuse 0.5
                $"<radialGradient id=\"item-{index}\" cx=\"%.3f{cx}\" cy=\"%.3f{cy}\" r=\"%.3f{r}\" fx=\"%.3f{cx}\" fy=\"%.3f{cy}\" gradientTransform=\"matrix(1,0,0,1,0,0)\" gradientUnits=\"userSpaceOnUse\"><stop offset=\"0.00\" stop-color=\"{startColor}\"/><stop offset=\"1.00\" stop-color=\"{stopColor}\"/></radialGradient>"
    
    /// <summary>
    /// Shape is a collection of supported shapes to draw in two-dimensional Euclidean space as SVG XML objects.
    /// </summary>
    type Shape =
        | Circle of int * Color * Circle * float 
        | Polygon of int * Color * Point2D list * float 
        | BondPolygon of int * Color * Point2D * Point2D * Point2D * Point2D * float // Bonds have a different gradient fill for Glossy style.
        | Line of int * Color * Point2D * Point2D * float 
        with
        member this.ToSvg (style : ArtStyle) =
            match this with 
            | Circle (index, color, circle, opacity) -> 
                let x = circle.Center.X
                let y = circle.Center.Y
                
                match style with
                | Cartoon ->
                    $"<circle cx=\"%.3f{x}\" cy=\"%.3f{y}\" r=\"%.3f{circle.Radius}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.05\" fill-opacity=\"{opacity}\" stroke-opacity=\"{opacity}\"/>"
                
                | Glossy ->
                    $"<circle class=\"item-{index}\" cx=\"%.3f{x}\" cy=\"%.3f{y}\" r=\"%.3f{circle.Radius}\" fill-opacity=\"{opacity}\"/>"
                
            | Polygon (index, color, points, opacity) ->
                let pointsStr = points |> List.map (fun p -> $"%.3f{p.X},%.3f{p.Y}") |> String.concat " "
                
                match style with
                | Cartoon -> 
                    $"<polygon points=\"{pointsStr}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.05\" stroke-linejoin=\"round\" fill-opacity=\"{opacity}\" stroke-opacity=\"{opacity}\"/>"
                    
                | Glossy ->
                    $"<polygon class=\"item-{index}\" points=\"{pointsStr}\" fill-opacity=\"{opacity}\"/>"
                    
            | BondPolygon (index, color, p1, p2, p3, p4, opacity) -> 
                let r1 = (p1.Dist p2) / 2.0
                let r2 = (p3.Dist p4) / 2.0
                
                let path = $"M %.3f{p1.X} %.3f{p1.Y} A %.3f{r1} %.3f{r1} 0 0 1 %.3f{p2.X} %.3f{p2.Y} L %.3f{p3.X} %.3f{p3.Y} A %.3f{r2} %.3f{r2} 0 0 1 %.3f{p4.X} %.3f{p4.Y} Z"
                
                match style with
                | Cartoon ->
                    $"<path d=\"{path}\" fill=\"{color}\" style=\"stroke:black;stroke-width:0.05\" stroke-linejoin=\"round\" fill-opacity=\"{opacity}\" stroke-opacity=\"{opacity}\"/>"
                
                | Glossy ->
                    $"<path class=\"item-{index}\" d=\"{path}\" fill-opacity=\"{opacity}\"/>"
            
            | Line (index, color, p1, p2, opacity) ->
                let x1 = p1.X
                let y1 = p1.Y
                let x2 = p2.X
                let y2 = p2.Y
                $"<line class=\"item-{index}\" x1=\"%.3f{x1}\" y1=\"%.3f{y1}\" x2=\"%.3f{x2}\" y2=\"%.3f{y2}\" stroke=\"{color}\" stroke-width=\"0.1\" stroke-linecap=\"round\" stroke-opacity=\"{opacity}\"/>"
        
        member this.Fill (style : ArtStyle) : Fill option =
            match style with
            | Cartoon -> None
            | Glossy ->
                match this with
                | Circle (index, color, circle, _) ->
                    (index, circle.Center, circle.Radius, color) |> RadialGradient |> Some 
                
                | Polygon (index, color, points, _) ->
                    let centroid = calcCentroid points
                    let maxDist = points |> List.map (fun p -> p.Dist centroid) |> List.max
                    (index, centroid, maxDist, color) |> RadialGradient |> Some    
                
                // | BondPolygon (index, color, points) ->    
                | BondPolygon (index, color, p1, p2, p3, p4, _) ->
                    let midPoint1 = p1.Midpoint p3
                    let midPoint2 = p2.Midpoint p4
                    (index, midPoint1, midPoint2, color) |> LinearGradient |> Some

                | _ -> None
        
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
    type SVG =
        { Header : Header
          ID : string
          ViewBox : ViewBox
          Objects : Shape list
          Style : ArtStyle }
        with
        override this.ToString () =            
            // Concatenate definitions, objects, and header strings. 
            this.Header.ToString() + this.Body() 
            
        member this.Body() =
            // Format style references.
            let formatStyleReference (o: Shape) : string =
                match o with
                | Circle (index, _, _, _)
                | Polygon (index, _, _, _)
                | BondPolygon (index, _, _, _, _, _, _)
                | Line (index, _, _, _, _) ->
                    $".item-{index}{{fill:url(#item-{index});}}"
            
            // Get style references.
            let styles = this.Objects |> List.map (fun x -> formatStyleReference x) |> String.concat "\n"
            
            // Convert all objects to a single string.
            let objs = this.Objects |> List.map (fun x -> x.ToSvg this.Style) |> String.concat "\n"
            
            // Combine items into SVG body.
            match this.Style with
            | Cartoon ->
                $"<svg id=\"{this.ID}\" xmlns=\"http://www.w3.org/2000/svg\" {this.ViewBox.ToString()}>\n{objs}\n</svg>"
                
            | Glossy ->
                let defs : string = this.Objects |> List.map (fun o -> o.Fill(this.Style).ToString()) |> String.concat "\n"
                 
                $"\n<svg id=\"{this.ID}\" xmlns=\"http://www.w3.org/2000/svg\" {this.ViewBox.ToString()}>\n<defs>\n<style>\n{styles}\n</style>\n{defs}\n</defs>\n{objs}\n</svg>"
            
module Drawing =
    
    open Svg
    
    /// <summary>
    /// Drawing options.
    /// </summary>
    type DrawingOptions =
        { ViewBox : ViewBox option
          ArtStyle : ArtStyle
          ModelStyle : ModelStyle
          DisplayHydrogenAtoms : bool 
          Resolution : int }
        with
        static member New () =
            { ViewBox = None
              ArtStyle = Cartoon 
              ModelStyle = SpaceFilling
              DisplayHydrogenAtoms = false
              Resolution = 40 }