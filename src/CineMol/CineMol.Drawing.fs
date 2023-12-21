module CineMol.Drawing

open CineMol.Types
open CineMol.Types.Chem
open Types.Geometry
open Types.Svg 
open Types.Drawing

/// <summary>
/// Driver code for rotating atoms in molecule.
/// </summary>
let rotate (mol: Molecule) (axis: Axis) (rad: float) =
    let rotateAtom (atom : Atom) = { atom with Position = atom.Position.Rotate axis rad }
    { mol with Atoms = List.map (fun atom -> rotateAtom atom) mol.Atoms }

/// <summary>
/// Algorithm for finding convex hull of a set of points in 2D.
/// </summary>
let quickHull2D (points : Point2D list) : int list =
    // Find index of point with minimum X coordinate.
    let minIndex : int =
        points
        |> List.mapi (fun i point -> (i, point))
        |> List.minBy (fun (_, point) -> point.X)
        |> fst
        
    // Find index of point with maximum X coordinate.
    let maxIndex : int =
        points
        |> List.mapi (fun i point -> (i, point))
        |> List.maxBy (fun (_, point) -> point.X)
        |> fst
        
    // Create a range.
    let createRange (length : int) : int list = [ 0 .. 1 .. length - 1 ]
        
    let rec processPoints (pointList : Point2D list, indexList : int list, a : int, b : int) : int list =
        let W : Vector2D = pointList[b].CreateVector pointList[a]
        
        let signedDist : (int * float) list =
            indexList
            |> List.map (fun i -> (i, pointList[i]))
            |> List.map (fun (i, point) ->
                let V : Vector2D = point.CreateVector pointList[a]
                let dist : float = V.Cross W
                (i, dist))
            |> List.filter (fun (i, dist) -> dist > 0.0 && i <> a && i <> b)
        
        match signedDist with
        | [] -> [ a; b ]
        | _ ->
            let maxDistIndex : int = signedDist |> List.maxBy (fun (_, dist) -> dist) |> fst
            let newIndexList : int list = List.map fst signedDist
            
            let left = processPoints (pointList, newIndexList, a, maxDistIndex)
            let right = processPoints (pointList, newIndexList, maxDistIndex, b)
            
            (left |> List.rev |> List.tail |> List.rev) @ right
    
    let left = processPoints (points, createRange points.Length, minIndex, maxIndex)
    let right = processPoints (points, createRange points.Length, maxIndex, minIndex)
   
    (left |> List.rev |> List.tail |> List.rev) @ (right |> List.rev |> List.tail |> List.rev) 
 
/// <summary>
/// Make a polygon from visible part of atom.
/// </summary>
let makePolygon (currAtom : Atom) (intersectsWith : Atom list) (resolution : int) : Point2D list =
    let sphere : Sphere = { Center = currAtom.Position; Radius = currAtom.Radius }
    
    let otherSpheres : Sphere list =
        intersectsWith |> List.map (fun atom -> { Center = atom.Position; Radius = atom.Radius })
        
    let notInsideOtherSphere (point : Point3D) : bool =
        otherSpheres |> List.forall (fun otherSphere -> otherSphere.Center.Dist point > otherSphere.Radius)
        
    // let quadrants : Point3D list list = sphere.PointsOnSphere resolution
    //
    // let points : Point2D list =
    //     quadrants
    //     |> List.map (fun quadrant -> quadrant |> List.filter notInsideOtherSphere)
    //     |> List.head
    //     |> List.map (fun point -> { X = point.X; Y = point.Y })
    
    let pointsOnSphere = sphere.PointsOnSphere resolution
    
    let points : Point2D list =
        pointsOnSphere 
        |> List.filter notInsideOtherSphere
        |> List.map (fun point -> { X = point.X; Y = point.Y })
        
    let indices : int list = quickHull2D points
    indices |> List.map (fun i -> points[i])
 
/// <summary>
/// Driver code for translating atoms in molecule.
/// </summary> 
let atomToSvg (prevAtoms : Atom list) (currAtom : Atom) (resolution : int) (opacity : float) : Shape =
    let currAtomGeom : Sphere = { Center = currAtom.Position; Radius = currAtom.Radius }
    
    let intersectsWith =
        prevAtoms
        |> List.filter (fun atom -> { Center = atom.Position; Radius = atom.Radius }.Intersects currAtomGeom)
        
    match intersectsWith with
    | [] ->
        let position : Point2D = { X = currAtom.Position.X; Y = currAtom.Position.Y }
        Shape.Circle (currAtom.Index, currAtom.Color, { Center = position; Radius = currAtom.Radius }, opacity )
        
    | _ ->
        let points : Point2D list = makePolygon currAtom intersectsWith resolution 
        Shape.Polygon (currAtom.Index, currAtom.Color, points, opacity)
 
/// <summary>
/// Driver code for translating bonds in molecule.
/// </summary>
let calculatePerpendicularPoints (p1 : Point3D) (p2 : Point3D) (radius : float) (numPoints : int) : Point3D list =
    let dx = p2.X - p1.X
    let dy = p2.Y - p1.Y
    let len = sqrt (dx ** 2.0 + dy ** 2.0)
    let ux = dx / len
    let uy = dy / len
    let step = 2.0 * radius / float (numPoints - 1)
    [
        for i in 0 .. numPoints - 1 ->
            let offset = float i * step - radius
            { X = p1.X + offset * uy; Y = p1.Y - offset * ux; Z = p1.Z }
    ]
    
let drawSingleBond (width : float) (index1 : int) (index2 : int) (color1 : Color) (color2 : Color) (p1 : Point3D) (p2 : Point3D) (opacity1 : float) (opacity2 : float) : Shape list =
    // Move p1 and p2 both a bit closer to each other, without using MoveTowards.
    let p1 = p1 + (p2 - p1).Mul(0.3)
    let p2 = p2 + (p1 - p2).Mul(0.3)
   
    let points1 = calculatePerpendicularPoints p1 p2 width 2
    let midPoints = calculatePerpendicularPoints (p1.Midpoint p2) p2 width 2 
    let points2 = calculatePerpendicularPoints p2 p1 width 2
    
    let b1, b2, b3, b4 : Point2D * Point2D * Point2D * Point2D =
        match points1 @ (midPoints |> List.rev) with
        | [ b1; b2; b3; b4 ] -> { X = b1.X; Y = b1.Y }, { X = b2.X; Y = b2.Y }, { X = b3.X; Y = b3.Y }, { X = b4.X; Y = b4.Y }
        | _ -> { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 } // TODO
    
    let b5, b6, b7, b8 : Point2D * Point2D * Point2D * Point2D =
        match midPoints @ points2 with
        | [ b5; b6; b7; b8 ] -> { X = b5.X; Y = b5.Y }, { X = b6.X; Y = b6.Y }, { X = b7.X; Y = b7.Y }, { X = b8.X; Y = b8.Y }
        | _ -> { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 }, { X = 0.0; Y = 0.0 } // TODO
    
    [ (index1, color1, b1, b2, b3, b4, opacity1) |> Shape.BondPolygon; (index2, color2, b5, b6, b7, b8, opacity2) |> Shape.BondPolygon ]

let bondToSvg (beginAtom : Atom) (endAtom : Atom) (bond : Bond) : Shape list =
    // Most negative Z is p1, most positive Z is p2, midpoint is pMid.
    let a1, a2 = if beginAtom.Position.Z > endAtom.Position.Z then endAtom, beginAtom else beginAtom, endAtom
    let p1 = a1.Position
    let p2 = a2.Position
    
    let index1, index2, opacity1, opacity2 =
        match bond.BeginAtomIndex, bond.EndAtomIndex with
        | x1, x2 when x1 = a1.Index && x2 = a2.Index -> bond.BeginIndex, bond.EndIndex, a1.Opacity, a2.Opacity
        | x1, x2 when x1 = a2.Index && x2 = a1.Index -> bond.EndIndex, bond.BeginIndex, a2.Opacity, a1.Opacity
        | _ -> 0, 0, 1.0, 1.0 // TODO: Catch this error properly.
    
    match bond.Type with
    | Single | Aromatic -> drawSingleBond 0.2 index1 index2 a1.Color a2.Color p1 p2 opacity1 opacity2
        
    | Double ->
        let radius = bond.Radius / 3.0
        let points1 = calculatePerpendicularPoints p1 p2 radius 2
        let points2 = calculatePerpendicularPoints p2 p1 radius 2 |> List.rev
        drawSingleBond 0.1 index1 index2 a1.Color a2.Color points1[0] points2[0] opacity1 opacity2 @
        drawSingleBond 0.1 index1 index2 a1.Color a2.Color points1[1] points2[1] opacity1 opacity2
        
    | Triple ->
        let radius = bond.Radius / 3.0
        let points1 = calculatePerpendicularPoints p1 p2 radius 3
        let points2 = calculatePerpendicularPoints p2 p1 radius 3 |> List.rev
        drawSingleBond 0.05 index1 index2 a1.Color a2.Color points1[0] points2[0] opacity1 opacity2 @
        drawSingleBond 0.05 index1 index2 a1.Color a2.Color points1[1] points2[1] opacity1 opacity2 @
        drawSingleBond 0.05 index1 index2 a1.Color a2.Color points1[2] points2[2] opacity1 opacity2
        
let bondToWire (beginAtom : Atom) (endAtom : Atom) : Shape list =
    // Most negative z is p1, most positive z is p2, midpoint is pMid.
    let a1, a2 = if beginAtom.Position.Z > endAtom.Position.Z then endAtom, beginAtom else beginAtom, endAtom
    let p1 = a1.Position
    let p2 = a2.Position
    let midPoint = p1.Midpoint p2
    
    let p1 : Point2D = { X = p1.X; Y = p1.Y }
    let p2 : Point2D = { X = p2.X; Y = p2.Y }
    let midPoint : Point2D = { X = midPoint.X; Y = midPoint.Y }
    
    // Index actually doesn't matter because the wire frame has a fixed style (since it's just lines).
    [ (a1.Index, a1.Color, p1, midPoint, a1.Opacity) |> Shape.Line; (a2.Index, a2.Color, midPoint, p2, a2.Opacity) |> Shape.Line ]
 
/// <summary>
/// Driver code for creating SVG for molecule.
/// </summary>
let draw (mol: Molecule) (options: DrawingOptions) : SVG * DrawingOptions =
    // Set origin.
    let origin = Axis.Origin()
    
    // Set view box.
    let viewBox  =
        let margin = 5.0
        
        let offset =
            mol.Atoms
            |> List.map (fun atom -> atom.Position.Dist origin)
            |> List.max
            |> (+) margin
            
        let viewBox =
            match options.ViewBox with
            | Some viewBox -> viewBox
            | None ->
                { MinX = -offset
                  MinY = -offset
                  Width = offset * 2.0
                  Height = offset * 2.0 }
                
        viewBox
        
    // Drawing style dictates if and how the objects are clipped and exactly drawn.
    let objs : Shape list =
        
        // Filter out hydrogen atoms.
        let atoms =
            match options.DisplayHydrogenAtoms with
            | true -> mol.Atoms
            | false -> mol.Atoms |> List.filter (fun atom -> atom.Type <> AtomType.H)
        
        match options.ModelStyle with        
        | SpaceFilling ->
            let atoms = atoms |> List.sortBy (fun atom -> atom.Position.Z) |> List.rev
    
            let rec processAtoms (prevAtoms : Atom list) (shapes : Shape list) : Shape list =
                match prevAtoms with
                | [] -> shapes
                | currAtom :: prevAtoms ->
                    let shape : Shape = atomToSvg prevAtoms currAtom options.Resolution currAtom.Opacity 
                    processAtoms prevAtoms (shape :: shapes)
                    
            match atoms with
            | [] -> []
            | atoms -> processAtoms atoms []
            
        | BallAndStick ->
            let atoms = atoms |> List.sortBy (fun atom -> atom.Position.Z) |> List.rev
            
            // BallAndStick has smaller atom radii.
            let atomRadiusCorrection = 0.3
            let atoms = atoms |> List.map (fun atom -> { atom with Radius = atom.Radius * atomRadiusCorrection })
            
            let rec processAtoms (prevAtoms : Atom list) (shapes : Shape list) : Shape list =
                match prevAtoms with
                | [] -> shapes
                | currAtom :: prevAtoms ->
                    let shape : Shape = atomToSvg prevAtoms currAtom options.Resolution currAtom.Opacity
                    
                    // Get bonds connected to atoms that have not yet been converted to an SVG shape.
                    let prevAtomIndices : int list = prevAtoms |> List.map (fun atom -> atom.Index)
                    let bondsToDraw =
                        mol.GetBonds options.DisplayHydrogenAtoms currAtom.Index
                        |> List.filter (fun bond ->
                            not (prevAtomIndices |> List.contains bond.BeginAtomIndex)
                            && not (prevAtomIndices |> List.contains bond.EndAtomIndex))
                        
                    // Convert bonds to SVG shapes.
                    let bondShapes : Shape list =
                        bondsToDraw
                        |> List.map (fun bond ->
                            let beginAtom = mol.GetAtom bond.BeginAtomIndex
                            let endAtom = mol.GetAtom bond.EndAtomIndex
                            match beginAtom, endAtom with
                            | Some b, Some e ->
                                let b = { b with Radius = b.Radius * atomRadiusCorrection }
                                let e = { e with Radius = e.Radius * atomRadiusCorrection }
                                bondToSvg b e bond  
                            | _ -> [])
                        |> List.concat
                    
                    processAtoms prevAtoms (shape :: bondShapes @ shapes)
            
            match atoms with
            | [] -> []
            | atoms -> processAtoms atoms []
            
        | WireFrame ->
            let atoms = atoms |> List.sortBy (fun atom -> atom.Position.Z) |> List.rev
            
            let rec processAtoms (prevAtoms : Atom list) (shapes : Shape list) : Shape list =
                match prevAtoms with
                | [] -> shapes
                | currAtom :: prevAtoms ->
                    
                    // Get bonds connected to atoms that have not yet been converted to an SVG shape.
                    let prevAtomIndices : int list = prevAtoms |> List.map (fun atom -> atom.Index)
                    let bondsToDraw =
                        mol.GetBonds options.DisplayHydrogenAtoms currAtom.Index
                        |> List.filter (fun bond ->
                            not (prevAtomIndices |> List.contains bond.BeginAtomIndex)
                            && not (prevAtomIndices |> List.contains bond.EndAtomIndex))
                        
                    // Convert bonds to SVG shapes.
                    let bondShapes : Shape list =
                        bondsToDraw
                        |> List.map (fun bond ->
                            let beginAtom = mol.GetAtom bond.BeginAtomIndex
                            let endAtom = mol.GetAtom bond.EndAtomIndex
                            match beginAtom, endAtom with
                            | Some b, Some e -> bondToWire b e  
                            | _ -> [])
                        |> List.concat
                    
                    processAtoms prevAtoms (bondShapes @ shapes)
            
            match atoms with
            | [] -> []
            | atoms -> processAtoms atoms []

    { Header = Header.New(); ID = "model"; ViewBox = viewBox; Objects = objs; Style = options.ArtStyle }, options           