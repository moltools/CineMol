module Client.CineMol.Styles


type RGB = int * int * int
type Alpha = float
type DiffusionGradient = Alpha[]

type Color = { RGB: RGB; Alpha: Alpha }
    with
    member x.Diffuse (alpha : float) : float * float * float =
        let R, G, B = x.RGB
        (float R) * alpha, (float G) * alpha, (float B * alpha)

let atomColorGradient: DiffusionGradient = [| 1.0; 0.89; 0.66; 0.34; 0.0 |]

type AtomType =
    | H  | He | Li | Be | B  | C  | N  | O  | F  | Ne | Na | Mg | Al | Si
    | P  | S  | Cl | Ar | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni
    | Cu | Zn | Ga | Ge | As | Se | Br | Kr | Rb | Sr | Y  | Zr | Nb | Mo
    | Tc | Ru | Rh | Pd | Ag | Cd | In | Sn | Sb | Te | I  | Xe | Cs | Ba
    | La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb
    | Lu | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg | Tl | Pb | Bi | Po
    | At | Rn | Fr | Ra | Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf
    | Es | Fm | Md | No | Lr | Rf | Db | Sg | Bh | Hs | Mt | Ds | Rg | Cn
    | Nh | Fl | Mc | Lv | Ts | Og | Unknown

type AtomColorStyle = | CPK

let getAtomColor (style: AtomColorStyle) (atomType: AtomType) : Color =
    match style with
    | CPK ->
        match atomType with
        | H  -> 255, 255, 255 | He -> 255, 255, 192 | Li -> 204, 128, 255
        | Be -> 194, 255, 000 | B  -> 255, 181, 181 | C  -> 144, 144, 144
        | N  -> 048, 080, 248 | O  -> 255, 013, 013 | F  -> 144, 224, 080
        | Ne -> 179, 227, 245 | Na -> 171, 092, 242 | Mg -> 138, 255, 000
        | Al -> 191, 166, 166 | Si -> 240, 200, 160 | P  -> 255, 128, 000
        | S  -> 255, 255, 048 | Cl -> 031, 240, 031 | Ar -> 128, 209, 227
        | K  -> 143, 064, 212 | Ca -> 061, 255, 000 | Sc -> 230, 230, 230
        | Ti -> 191, 194, 199 | V  -> 166, 166, 171 | Cr -> 138, 153, 199
        | Mn -> 156, 122, 199 | Fe -> 224, 102, 051 | Co -> 240, 144, 160
        | Ni -> 080, 208, 080 | Cu -> 200, 128, 051 | Zn -> 125, 128, 176
        | Ga -> 194, 143, 143 | Ge -> 102, 143, 143 | As -> 189, 128, 227
        | Se -> 255, 161, 000 | Br -> 166, 041, 041 | Kr -> 092, 184, 209
        | Rb -> 112, 046, 176 | Sr -> 000, 255, 000 | Y  -> 148, 255, 255
        | Zr -> 148, 224, 224 | Nb -> 115, 194, 201 | Mo -> 084, 181, 181
        | Tc -> 059, 158, 158 | Ru -> 036, 143, 143 | Rh -> 010, 125, 140
        | Pd -> 000, 105, 133 | Ag -> 192, 192, 192 | Cd -> 255, 217, 143
        | In -> 166, 117, 115 | Sn -> 102, 128, 128 | Sb -> 158, 099, 181
        | Te -> 212, 122, 000 | I  -> 148, 000, 148 | Xe -> 066, 158, 176
        | Cs -> 087, 023, 143 | Ba -> 000, 201, 000 | La -> 112, 212, 255
        | Ce -> 255, 255, 199 | Pr -> 217, 255, 199 | Nd -> 199, 255, 199
        | Pm -> 163, 255, 199 | Sm -> 143, 255, 199 | Eu -> 097, 255, 199
        | Gd -> 069, 255, 199 | Tb -> 048, 255, 199 | Dy -> 031, 255, 199
        | Ho -> 000, 255, 156 | Er -> 000, 230, 117 | Tm -> 000, 212, 082
        | Yb -> 000, 191, 056 | Lu -> 000, 171, 036 | Hf -> 077, 194, 255
        | Ta -> 077, 166, 255 | W  -> 033, 148, 214 | Re -> 038, 125, 171
        | Os -> 038, 102, 150 | Ir -> 023, 084, 135 | Pt -> 208, 208, 224
        | Au -> 255, 209, 035 | Hg -> 184, 184, 208 | Tl -> 166, 084, 077
        | Pb -> 087, 089, 097 | Bi -> 158, 079, 181 | Po -> 171, 092, 000
        | At -> 117, 079, 069 | Rn -> 066, 130, 150 | Fr -> 066, 000, 102
        | Ra -> 000, 125, 000 | Ac -> 112, 171, 250 | Th -> 000, 186, 255
        | Pa -> 000, 161, 255 | U  -> 000, 143, 255 | Np -> 000, 128, 255
        | Pu -> 000, 107, 255 | Am -> 084, 092, 242 | Cm -> 120, 092, 227
        | Bk -> 138, 079, 227 | Cf -> 161, 054, 212 | Es -> 179, 031, 212
        | Fm -> 179, 031, 186 | Md -> 179, 013, 166 | No -> 189, 013, 135
        | Lr -> 199, 000, 102 | Rf -> 204, 000, 089 | Db -> 209, 000, 079
        | Sg -> 217, 000, 069 | Bh -> 224, 000, 056 | Hs -> 230, 000, 046
        | Mt -> 235, 000, 038 | _  -> 000, 000, 000
        |> (fun rgb -> { RGB = rgb; Alpha = 1.0 })

type AtomGeomStyle = | Default
type Radius = float

let getAtomRadius (style: AtomGeomStyle) (atomType: AtomType) : Radius =
    match style with
    | Default ->
        // Atom radii in pico meters (based on quantum mechanical wave functions).
        match atomType with
        | H  -> 037 | Li -> 152 | Be -> 112 | B  -> 088 | C  -> 077
        | N  -> 070 | O  -> 066 | F  -> 064 | Na -> 186 | Mg -> 160
        | Al -> 143 | Si -> 117 | P  -> 110 | S  -> 104 | Cl -> 099
        | K  -> 231 | Ca -> 197 | Sc -> 160 | Ti -> 146 | V  -> 131
        | Cr -> 125 | Mn -> 129 | Fe -> 126 | Co -> 126 | Ni -> 124
        | Cu -> 128 | Zn -> 133 | Ga -> 122 | Ge -> 122 | As -> 121
        | Se -> 117 | Br -> 114 | Rb -> 241 | Sr -> 215 | Y  -> 180
        | Zr -> 157 | Nb -> 143 | Mo -> 136 | Tc -> 130 | Ru -> 133
        | Rh -> 134 | Pd -> 138 | Ag -> 144 | Cd -> 149 | In -> 162
        | Sn -> 140 | Sb -> 141 | Te -> 137 | I  -> 133 | Cs -> 262
        | Ba -> 217 | Hf -> 157 | Ta -> 143 | W  -> 137 | Re -> 137
        | Os -> 134 | Ir -> 135 | Pt -> 138 | Au -> 144 | Hg -> 150
        | Tl -> 171 | Pb -> 175 | Bi -> 146 | Po -> 140 | At -> 140
        | _  -> 077
        |> float

let normalizeRadius (style: AtomGeomStyle) (radius: Radius) : Radius =
    match style with
    // Default atom radius is normalized by dividing by hydrogen atom radius.
    | Default -> radius / (getAtomRadius Default AtomType.C)