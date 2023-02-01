module CineMol.Drawing

open System 

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
        let objs: Shape list =
            adjAtoms
            |> List.map (fun (Atom2D ({ Index = index; Type = _; Color = color }, c, Radius r)) ->
                (index, color, Circle2D(c, Radius (r / 5.0))) |> Circle)
        
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
        
    | Tube ->
        // Tube model depiction.
        // TODO: add specular, and add masks to tips of bonds to make them look round
        // TODO: make sure SVG elements resize with chaning perspective 
        
        // Create atom atom index to atom look-up map. 
        let getAtom = adjAtoms |> List.map (fun a -> (a.GetInfo().Index, a)) |> Collections.Map 
        
        // Find bonds connected to atom with certain atom index.
        let getBonds idx = mol.Bonds |> List.filter (fun (Bond b) -> b.BeginAtomIndex = idx)
        
        // Keep track which atoms have been drawn (as joints between wires) to prevent drawing the same atom multiple times.
        let mutable drawn = []
        
        // Thickness of wire-frame model in Angstrom.
        let width = 0.1

        // Draw bonds as tubes.
        let drawBonds (Atom2D(bInfo, bCenter, _)) (Atom2D(eInfo, eCenter, _)) (Bond bond) width =
            let slopePerpendicular = bCenter.SlopePerpendicular eCenter
            
            let translation width =
                let t = width / Math.Sqrt (1.0 + Math.Pow(slopePerpendicular, 2.0))
                { X = t; Y = slopePerpendicular * t }
                
            let elem idx color bL bR eR eL = (idx, color, Types.Geometry.Quadrangle (bL, bR, eR, eL)) |> Quadrangle
            
            let drawTube (b: Point2D) bIdx bCol e eIdx eCol adj =
                let m = b.Midpoint e 
                [ elem bIdx bCol (b + adj) (b - adj) (m - adj) (m + adj)
                  elem eIdx eCol (e + adj) (e - adj) (m - adj) (m + adj) ]
                
            match bond.Type with
            | Single | Aromatic -> drawTube bCenter bInfo.Index bInfo.Color eCenter eInfo.Index eInfo.Color (translation width)
            | Double ->
                let sep, newWidth = width / 2.0, width / 3.0
                let adj = translation sep
                let bL, bR, eL, eR = bCenter + adj, bCenter - adj, eCenter + adj, eCenter - adj
                drawTube bL bInfo.Index bInfo.Color eL eInfo.Index eInfo.Color (translation newWidth) @
                drawTube bR bInfo.Index bInfo.Color eR eInfo.Index eInfo.Color (translation newWidth)
            | Triple ->
                let sep, newWidth = width / 2.0, width / 5.0
                let adj = translation sep
                let bL, bR, eL, eR = bCenter + adj, bCenter - adj, eCenter + adj, eCenter - adj
                drawTube bL bInfo.Index bInfo.Color eL eInfo.Index eInfo.Color (translation newWidth) @
                drawTube bCenter bInfo.Index bInfo.Color eCenter eInfo.Index eInfo.Color (translation newWidth) @
                drawTube bR bInfo.Index bInfo.Color eR eInfo.Index eInfo.Color (translation newWidth) 
                        
        // Draw objects.
        let objs = [
            for bAtom in adjAtoms do
                yield Circle (bAtom.GetInfo().Index, bAtom.GetInfo().Color, Circle2D (bAtom.GetCenter(), Radius width))
                drawn <- drawn @ [ bAtom.GetInfo().Index ]
                for bond in getBonds (bAtom.GetInfo().Index) do
                    if not (List.contains (bond.Unwrap().EndAtomIndex) drawn) then
                        match getAtom.TryFind (bond.Unwrap().EndAtomIndex) with
                        | None -> () | Some eAtom -> for bondTube in drawBonds eAtom bAtom bond width do yield bondTube ]
        
        { Header = Header.New(); ID = "Tube"; ViewBox = viewBox; Objects = objs }, options