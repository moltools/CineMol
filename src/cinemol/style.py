from abc import ABC, abstractmethod
from dataclasses import dataclass

class Color:
    def __init__(self, r: int, g: int, b: int) -> None: 
        self.r = max(0, min(255, r))
        self.g = max(0, min(255, g))
        self.b = max(0, min(255, b))

    def __str__(self) -> str:
        return f"rgb({self.r},{self.g},{self.b}"
    
    def to_hex(self) -> str:
        return f"#{self.r:02x}{self.g:02x}{self.b:02x}"
    
    def diffuse(self, alpha: float) -> "Color":
        return Color(self.r * alpha, self.g * alpha, self.b * alpha)
    
class AtomColoringScheme(ABC):
    @abstractmethod
    def get_color(self, atom_symbol: str) -> Color:
        ...

class AtomRadiusScheme(ABC):
    @abstractmethod
    def to_angstrum(self, atom_symbol: str) -> float:
        ...

class CoreyPaulingKoltungAtomColor(AtomColoringScheme):
    """
    Corey-Pauling-Koltun coloring convention for atoms. 
    
    Source: https://en.wikipedia.org/wiki/CPK_coloring
    """
    H  = Color(255, 255, 255) # White
    C  = Color( 48,  48,  48) # Not quite black (need non-zero values to able to diffuse)
    N  = Color(  0,   0, 255) # Blue
    O  = Color(255,   0,   0) # Red
    P  = Color(255, 165,   0) # Orange
    S  = Color(255, 255,   0) # Yellow
    B  = Color(245, 245, 220) # Beige
    Br = Color(139,   0,   0) # Dark red
    I  = Color(148,   0, 211) # Dark violet
    Ti = Color(128, 128, 128) # Grey
    Fe = Color(255, 140,   0) # Dark orange
    F  = Color(  0, 128,   0) # Green
    Cl = Color(  0, 128,   0) # Green
    He = Color(  0, 255, 255) # Cyan
    Ne = Color(  0, 255, 255) # Cyan
    Ar = Color(  0, 255, 255) # Cyan
    Kr = Color(  0, 255, 255) # Cyan
    Xe = Color(  0, 255, 255) # Cyan
    Li = Color(238, 130, 238) # Violet
    Na = Color(238, 130, 238) # Violet
    K  = Color(238, 130, 238) # Violet
    Rb = Color(238, 130, 238) # Violet
    Cs = Color(238, 130, 238) # Violet
    Fr = Color(238, 130, 238) # Violet
    Be = Color(  0, 100,   0) # Dark green
    Mg = Color(  0, 100,   0) # Dark green
    Ca = Color(  0, 100,   0) # Dark green
    Sr = Color(  0, 100,   0) # Dark green
    Ba = Color(  0, 100,   0) # Dark green
    Ra = Color(  0, 100,   0) # Dark green

    def get_color(self, atom_symbol: str) -> Color:
        default = Color(255, 192, 203) # Pink  
        return getattr(self, atom_symbol, default)

class PubChemAtomRadius(AtomRadiusScheme):
    """
    Atomic radii (van der Waals) in picometer from PubChem.
    
    Source: https://pubchem.ncbi.nlm.nih.gov/periodic-table/#property=AtomicRadius 
    """
    H  = 120.0; He = 140.0; Li = 182.0; Be = 153.0; B  = 192.0; C  = 170.0
    N  = 155.0; O  = 152.0; F  = 135.0; Ne = 154.0; Na = 227.0; Mg = 173.0
    Al = 184.0; Si = 210.0; P  = 180.0; S  = 180.0; Cl = 175.0; Ar = 188.0
    K  = 275.0; Ca = 231.0; Sc = 211.0; Ti = 187.0; V  = 179.0; Cr = 189.0
    Mn = 197.0; Fe = 194.0; Co = 192.0; Ni = 163.0; Cu = 140.0; Zn = 139.0
    Ga = 187.0; Ge = 211.0; As = 185.0; Se = 190.0; Br = 183.0; Kr = 202.0
    Rb = 303.0; Sr = 249.0; Y  = 219.0; Zr = 186.0; Nb = 207.0; Mo = 209.0
    Tc = 209.0; Ru = 207.0; Rh = 195.0; Pd = 202.0; Ag = 172.0; Cd = 158.0
    In = 193.0; Sn = 217.0; Sb = 206.0; Te = 206.0; I  = 198.0; Xe = 216.0
    Cs = 343.0; Ba = 268.0; Lu = 221.0; Hf = 212.0; Ta = 217.0; W  = 210.0
    Re = 217.0; Os = 216.0; Ir = 202.0; Pt = 209.0; Au = 166.0; Hg = 209.0
    Tl = 196.0; Pb = 202.0; Bi = 207.0; Po = 197.0; At = 202.0; Rn = 220.0
    Fr = 348.0; Ra = 283.0

    def to_angstrum(self, atom_symbol: str) -> float:
        default = 170.0
        return getattr(self, atom_symbol, default) / 100

@dataclass
class Style:
    reference: str 
    color: Color

    def to_svg(self) -> str:
        return f".{self.reference}{{fill:{self.color.to_hex()}}}"