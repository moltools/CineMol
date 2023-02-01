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
        
        // Create atom atom index to atom look-up map. 
        let lookUp =
            adjAtoms
            |> List.map (fun atom ->
                let (Atom2D({ Index = Index idx; Type = _; Color = _ }, _, _)) = atom
                (idx, atom))
            |> Collections.Map 
        
        // Construct wire-fram from bonds.
        let objs: Shape list =
            mol.Bonds
            
            // Retrieve atoms which are connected by bonds to draw.
            |> List.map (fun bond ->
                let (Bond info) = bond 
                let (Index bIdx) = info.BeginAtomIndex
                let (Index eIdx) = info.EndAtomIndex
                match lookUp.TryFind bIdx, lookUp.TryFind eIdx with 
                | Some b, Some e -> Some (b, e, bond)
                | _ -> None)
            
            // Ignore bonds that refer to non-existing atom indices.
            |> List.choose id
            
            // Sort bonds based on begin atom from furthest away to nearest.
            |> List.sortBy (fun (Atom2D(_, c, _), _, _) -> -(c.Dist (pov.ToPoint2D())))
            
            // Draw bonds.
            |> List.map (fun (Atom2D(bInfo, bCenter, Radius bRadius), Atom2D(eInfo, eCenter, Radius eRadius), bond) ->
                // TODO: draw bonds
                // TODO: draw double, triple, and aromatic bonds
                // TODO: if atom is starting point for multiple bonds it is now drawn multiple times 
                (bInfo.Index, bInfo.Color, Circle2D (bCenter, Radius (bRadius / 10.0))) |> Circle)
        
        printf "Dm: %A" <| objs 
        
        { Header = Header.New(); ID = "WireFrame"; ViewBox = viewBox; Objects = objs }, options