module CineMol.Drawing

open Types.Fundamentals 
open Types.Geometry
open Types.Chem
open Types.Svg
open Types.Drawing
open Projection

/// <summary>
/// Driver code for rotating atoms in molecule.
/// </summary>
let rotate (mol: Molecule) (axis: Axis) (rad: float) =
    let rotateAtom (Atom3D (i, c, r): Atom3D) = Atom3D (i, c.Rotate axis rad , r)
    { mol with Atoms = List.map (fun atom -> rotateAtom atom) mol.Atoms }
    
/// <summary>
/// Driver code for creating SVG for molecule.
/// </summary>
let draw (mol: Molecule) (options: DrawingOptions) =
    // Set origin.
    let origin = Axis.Origin()
    
    // Set view box.
    let offset, viewBox  =
        let marginRatio = 2.0
        
        let offset =
            mol.Atoms
            |> List.map (fun (Atom3D (_, c, _)) -> c.Dist origin)
            |> List.max
            |> (*) marginRatio
            
        let viewBox =
            match options.ViewBox with
            | Some viewBox -> viewBox
            | None ->
                { MinX = -offset
                  MinY = -offset
                  Width = offset * marginRatio
                  Height = offset * marginRatio }
                
        offset, viewBox
        
    // Set point of view. The view box looks along the Z axis. 
    let pov = { X = 1E-5; Y = 1E-5; Z = offset }
    
    // Sort atoms based on distance atoms to point of view.
    let adjAtoms = mol.Atoms |> List.sortBy (fun (Atom3D (_, c, _)) -> -(c.Dist pov))
    
    // Reset atom radii based based on distance atom to point of view.
    let adjAtoms =
        adjAtoms 
        |> List.map (fun (Atom3D (i, c, Radius r)) ->
            Atom3D(i, c, Radius <| r * ((pov.Dist origin) / (pov.Dist c))))
        
    let adjAtoms =
        let project (p: Point3D) = project (Camera.New pov) pov offset p
        adjAtoms |> List.map (fun (Atom3D (i, c, r)) -> Atom2D (i, project c, r))
        
    // Drawing style dictates if and how the objects are clipped and exactly drawn.
    match options.Style with
    
    | BallAndStick ->
        // Ball-and-stick model depiction.
        // TODO: parse element 
        // TODO: add specular 
        
        // Convert elements to SVG objects.
        let objs: Shape list = []
        
        { Header = Header.New(); ID = "BallAndStick"; ViewBox = viewBox; Objects = objs }, options
        
    | SpaceFilling ->
        // Space-filling model depiction.
        // TODO: parse elements (incl. clipping)
        // TODO: add specular 
    
        // Convert elements to SVG objects.
        let objs: Shape list =
            adjAtoms
            |> List.map (fun (Atom2D ({ Index = index; Type = _; Color = color }, c, r)) ->
                (index, color, Circle2D(c, r)) |> Circle)
        
        { Header = Header.New(); ID = "SpaceFilling"; ViewBox = viewBox; Objects = objs }, options
        
    | WireFrame ->
        // Wire-frame model depiction.
        // TODO: parse elements
        // TODO: add specular 
        
        // Convert elements to SVG objects.
        let objs: Shape list = []
        
        { Header = Header.New(); ID = "WireFrame"; ViewBox = viewBox; Objects = objs }, options