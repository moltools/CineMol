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
///
/// </summary> 
let toSvg (prevAtoms : Atom list) (currAtom : Atom) : Shape =
    let position = { X = currAtom.Position.X; Y = currAtom.Position.Y }
    Shape.Circle (currAtom.Index, currAtom.Color, { Center = position; Radius = currAtom.Radius } )
 
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
            let atoms = mol.Atoms |> List.sortBy (fun atom -> atom.Position.Z) |> List.rev
    
            let rec processAtoms (prevAtoms : Atom list) (shapes : Shape list) : Shape list =
                match prevAtoms with
                | [] -> shapes
                | currAtom :: prevAtoms ->
                    let shape : Shape = toSvg prevAtoms currAtom
                    processAtoms prevAtoms (shape :: shapes)
                    
            match atoms with
            | [] -> []
            | atoms -> processAtoms atoms [] |> List.rev

    // { Header = Header.New(); ID = "model"; ViewBox = viewBox; Objects = objs; Masks = masks }, options
    { Header = Header.New(); ID = "model"; ViewBox = viewBox; Objects = objs }, options           