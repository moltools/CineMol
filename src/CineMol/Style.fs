module CineMol.Style

open CineMol.Types.Fundamentals
open CineMol.Types.Style 
open CineMol.Types.Chem

type AtomColorStyle =
    
    /// Corey-Pauling-Koltun coloring convention for atoms.
    /// Source: https://en.wikipedia.org/wiki/CPK_coloring  
    | CPK
    
    with
    member this.Color (atomType: AtomType) =
        match atomType with
        | H                           -> 255, 255, 255 // White
        | C                           ->   0,   0,   0 // Black 
        | N                           ->   0,   0, 255 // Blue
        | O                           -> 255,   0,   0 // Red
        | P                           -> 255, 165,   0 // Orange
        | S                           -> 255, 255,   0 // Yellow
        | B                           -> 245, 245, 220 // Beige
        | Br                          -> 139,   0,   0 // Dark red
        | I                           -> 148,   0, 211 // Dark violet
        | Ti                          -> 128, 128, 128 // Grey
        | Fe                          -> 255, 140,   0 // Dark orange
        | F  | Cl                     ->   0, 128,   0 // Green
        | He | Ne | Ar | Kr | Xe      ->   0, 255, 255 // Cyan
        | Li | Na | K  | Rb | Cs | Fr -> 238, 130, 238 // Violet
        | Be | Mg | Ca | Sr | Ba | Ra ->   0, 100,   0 // Dark green
        | _                           -> 255, 192, 203 // Pink 
        |> Color
        
type AtomRadius =
    
    /// Atomic radii (van der Waals) in pm from PubChem.
    /// Source: https://pubchem.ncbi.nlm.nih.gov/periodic-table/#property=AtomicRadius 
    | PubChem
    
    with
    member this.Radius (atomType: AtomType) =
        match atomType with
        | H  -> 120.0 | He -> 140.0 | Li -> 182.0 | Be -> 153.0 | B  -> 192.0 | C  -> 170.0
        | N  -> 155.0 | O  -> 152.0 | F  -> 135.0 | Ne -> 154.0 | Na -> 227.0 | Mg -> 173.0
        | Al -> 184.0 | Si -> 210.0 | P  -> 180.0 | S  -> 180.0 | Cl -> 175.0 | Ar -> 188.0
        | K  -> 275.0 | Ca -> 231.0 | Sc -> 211.0 | Ti -> 187.0 | V  -> 179.0 | Cr -> 189.0
        | Mn -> 197.0 | Fe -> 194.0 | Co -> 192.0 | Ni -> 163.0 | Cu -> 140.0 | Zn -> 139.0
        | Ga -> 187.0 | Ge -> 211.0 | As -> 185.0 | Se -> 190.0 | Br -> 183.0 | Kr -> 202.0
        | Rb -> 303.0 | Sr -> 249.0 | Y  -> 219.0 | Zr -> 186.0 | Nb -> 207.0 | Mo -> 209.0
        | Tc -> 209.0 | Ru -> 207.0 | Rh -> 195.0 | Pd -> 202.0 | Ag -> 172.0 | Cd -> 158.0
        | In -> 193.0 | Sn -> 217.0 | Sb -> 206.0 | Te -> 206.0 | I  -> 198.0 | Xe -> 216.0
        | Cs -> 343.0 | Ba -> 268.0 | Lu -> 221.0 | Hf -> 212.0 | Ta -> 217.0 | W  -> 210.0
        | Re -> 217.0 | Os -> 216.0 | Ir -> 202.0 | Pt -> 209.0 | Au -> 166.0 | Hg -> 209.0
        | Tl -> 196.0 | Pb -> 202.0 | Bi -> 207.0 | Po -> 197.0 | At -> 202.0 | Rn -> 220.0
        | Fr -> 348.0 | Ra -> 283.0
        |> Radius
        
// TODO: how to normalize atom radii im pm to viewer? Perhaps based on average single CC bond length?
        