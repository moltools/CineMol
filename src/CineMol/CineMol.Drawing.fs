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
        
        // Create atom atom index to atom look-up map. 
        let getAtom = adjAtoms |> List.map (fun a -> (a.GetInfo().Index, a)) |> Collections.Map 
        
        // Find bonds connected to atom with certain atom index.
        let getBonds idx = mol.Bonds |> List.filter (fun (Bond b) -> b.BeginAtomIndex = idx)
        
        // Keep track which atoms have been drawn (as joints between wires) to prevent drawing the same atom multiple times.
        let mutable drawn = []
        
        // Resize atoms to ratio for wire-frame model in Angstrom.
        let factor = 10.0

        // Draw bonds as tubes.
        let drawBonds (Atom2D(bInfo, bCenter, bRad)) (Atom2D(eInfo, eCenter, eRad)) (Bond bond) =
            let slopePerpendicular = bCenter.SlopePerpendicular eCenter
            
            let translation width =
                let t = width / Math.Sqrt (1.0 + Math.Pow(slopePerpendicular, 2.0))
                { X = t; Y = slopePerpendicular * t }
                
            let elem idx color bL bR eR eL = (idx, color, Types.Geometry.Quadrangle (bL, bR, eR, eL)) |> Quadrangle
            
            let drawTube (b: Point2D) bIdx bCol bWidth e eIdx eCol eWidth =
                let m = b.Midpoint e
                let bAdj, eAdj, mAdj = translation bWidth, translation eWidth, translation ((bWidth + eWidth) / 2.0)
                [ elem bIdx bCol (b + bAdj) (b - bAdj) (m - mAdj) (m + mAdj)
                  elem eIdx eCol (e + eAdj) (e - eAdj) (m - mAdj) (m + mAdj) ]              
                
            match bond.Type with
            | Single | Aromatic ->
                drawTube bCenter bInfo.Index bInfo.Color (bRad.Unwrap() / factor) eCenter eInfo.Index eInfo.Color (eRad.Unwrap() / factor)
            | Double ->
                let bSepAdj, bNewWidth = translation ((bRad.Unwrap() / factor) / 2.0), (bRad.Unwrap() / factor) / 3.0
                let eSepAdj, eNewWidth = translation ((eRad.Unwrap() / factor) / 2.0), (eRad.Unwrap() / factor) / 3.0
                let bL, eL = bCenter + bSepAdj, eCenter + eSepAdj
                let bR, eR = bCenter - bSepAdj, eCenter - eSepAdj
                drawTube bL bInfo.Index bInfo.Color bNewWidth eL eInfo.Index eInfo.Color eNewWidth @
                drawTube bR bInfo.Index bInfo.Color bNewWidth eR eInfo.Index eInfo.Color eNewWidth
            | Triple ->
                let bSepAdj, bNewWidth = translation ((bRad.Unwrap() / factor) / 2.0), (bRad.Unwrap() / factor) / 5.0
                let eSepAdj, eNewWidth = translation ((eRad.Unwrap() / factor) / 2.0), (eRad.Unwrap() / factor) / 5.0
                let bL, eL = bCenter + bSepAdj, eCenter + eSepAdj
                let bR, eR = bCenter - bSepAdj, eCenter - eSepAdj
                drawTube bL bInfo.Index bInfo.Color bNewWidth eL eInfo.Index eInfo.Color eNewWidth @
                drawTube bCenter bInfo.Index bInfo.Color bNewWidth eCenter eInfo.Index eInfo.Color eNewWidth @
                drawTube bR bInfo.Index bInfo.Color bNewWidth eR eInfo.Index eInfo.Color eNewWidth
            
        // Draw objects.
        let objs = [
            for bAtom in adjAtoms do
                yield Circle (bAtom.GetInfo().Index, bAtom.GetInfo().Color, Circle2D (bAtom.GetCenter(), Radius (bAtom.GetRadius().Unwrap() / factor)))
                drawn <- drawn @ [ bAtom.GetInfo().Index ]
                for bond in getBonds (bAtom.GetInfo().Index) do
                    if not (List.contains (bond.Unwrap().EndAtomIndex) drawn) then
                        match getAtom.TryFind (bond.Unwrap().EndAtomIndex) with
                        | None -> () | Some eAtom -> for bondTube in drawBonds eAtom bAtom bond do yield bondTube ]
        
        { Header = Header.New(); ID = "Tube"; ViewBox = viewBox; Objects = objs }, options