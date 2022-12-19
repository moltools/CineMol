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
        | H  -> 053 | He -> 031 | Li -> 167 | Be -> 112 | B  -> 087 | C  -> 067
        | N  -> 056 | O  -> 048 | F  -> 042 | Ne -> 038 | Na -> 190 | Mg -> 145
        | Al -> 118 | Si -> 111 | P  -> 098 | S  -> 088 | Cl -> 079 | Ar -> 071
        | K  -> 243 | Ca -> 194 | Sc -> 184 | Ti -> 176 | V  -> 171 | Cr -> 166
        | Mn -> 161 | Fe -> 156 | Co -> 152 | Ni -> 149 | Cu -> 145 | Zn -> 142
        | Ga -> 136 | Ge -> 125 | As -> 114 | Se -> 103 | Br -> 094 | Kr -> 088
        | Rb -> 265 | Sr -> 219 | Y  -> 212 | Zr -> 206 | Nb -> 198 | Mo -> 190
        | Tc -> 183 | Ru -> 178 | Rh -> 173 | Pd -> 169 | Ag -> 165 | Cd -> 155
        | In -> 156 | Sn -> 145 | Sb -> 133 | Te -> 123 | I  -> 115 | Xe -> 108
        | Cs -> 298 | Ba -> 253 | Hf -> 208 | Ta -> 200 | W  -> 193 | Re -> 188
        | Os -> 185 | Ir -> 180 | Pt -> 177 | Au -> 174 | Hg -> 171 | Tl -> 156
        | Pb -> 154 | Bi -> 143 | Po -> 135 | At -> 127 | Rn -> 120 | _  -> 031

let normalizeRadius (style: AtomGeomStyle) (radius: Radius) : Radius =
    match style with
    // Default atom radius is normalized by dividing by hydrogen atom radius.
    | Default -> radius / (getAtomRadius Default AtomType.H)