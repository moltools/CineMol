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
            let maxDistIndex : int = signedDist |> List.maxBy (fun (i, dist) -> dist) |> fst
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
let toSvg (prevAtoms : Atom list) (currAtom : Atom) (resolution : int) : Shape =
    let currAtomGeom : Sphere = { Center = currAtom.Position; Radius = currAtom.Radius }
    
    let intersectsWith =
        prevAtoms
        |> List.filter (fun atom -> { Center = atom.Position; Radius = atom.Radius }.Intersects currAtomGeom)
        
    match intersectsWith with
    | [] ->
        let position : Point2D = { X = currAtom.Position.X; Y = currAtom.Position.Y }
        Shape.Circle (currAtom.Index, currAtom.Color, { Center = position; Radius = currAtom.Radius } )
        
    | _ ->
        let points : Point2D list = makePolygon currAtom intersectsWith resolution 
        Shape.Polygon (currAtom.Index, currAtom.Color, points)
 
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
        match options.Style with        
        | SpaceFilling ->
            // Sort atoms by position on Z axis.
            let atoms = mol.Atoms |> List.sortBy (fun atom -> atom.Position.Z) |> List.rev
    
            let rec processAtoms (prevAtoms : Atom list) (shapes : Shape list) : Shape list =
                match prevAtoms with
                | [] -> shapes
                | currAtom :: prevAtoms ->
                    let shape : Shape = toSvg prevAtoms currAtom options.Resolution 
                    processAtoms prevAtoms (shape :: shapes)
                    
            match atoms with
            | [] -> []
            | atoms -> processAtoms atoms [] 

    // { Header = Header.New(); ID = "model"; ViewBox = viewBox; Objects = objs; Masks = masks }, options
    { Header = Header.New(); ID = "model"; ViewBox = viewBox; Objects = objs }, options           