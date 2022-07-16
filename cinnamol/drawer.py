from asyncio import transports
import numpy as np
from enum import auto, Enum
from dataclasses import dataclass
from typing import List, Tuple

from pyvista import Plotter, Sphere, Tube

from error_handling import error_handler


# =============================================================================
# General definitions and functionalities.
# =============================================================================
@dataclass 
class Point:
    """
    Defines point in 3D space.
    """
    x: float 
    y: float 
    z: float 

    def as_array(self) -> np.array:
        """
        Returns point as [x, y, z] array.
        """
        return np.array([self.x, self.y, self.z])


@dataclass
class Radius:
    """
    Defines the radius of a sphere.
    """
    radius: float


@dataclass 
class Color:
    """
    Defines a color in RGB format.
    """
    r: float
    g: float 
    b: float


# =============================================================================
# Atom definitions and functionalities.
# =============================================================================
class Symbol(Enum):
    """
    Allowed atom types.
    """
    H  = auto() # Hydrogen
    He = auto() # Helium
    Li = auto() # Lithium
    Be = auto() # Beryllium
    B  = auto() # Boron
    C  = auto() # Carbon 
    N  = auto() # Nitrogen
    O  = auto() # Oxygen
    F  = auto() # Fluorine
    Ne = auto() # Neon 
    Na = auto() # Sodium
    Mg = auto() # Magnesium
    Al = auto() # Aluminum
    Si = auto() # Silicon 
    P  = auto() # Phosphorus
    S  = auto() # Sulfur
    Cl = auto() # Chlorine 
    Ar = auto() # Argon
    K  = auto() # Potassium
    Ca = auto() # Calcium
    Sc = auto() # Scandium
    Ti = auto() # Titanium
    V  = auto() # Vanadium
    Cr = auto() # Chromium
    Mn = auto() # Manganese
    Fe = auto() # Iron
    Co = auto() # Cobalt
    Ni = auto() # Nickel
    Cu = auto() # Copper
    Zn = auto() # Zinc
    Ga = auto() # Gallium
    Ge = auto() # Germanium
    As = auto() # Arsenic
    Se = auto() # Selenium
    Br = auto() # Bromine
    Kr = auto() # Krypton
    Rb = auto() # Rubidium
    Sr = auto() # Strontium
    Y  = auto() # Yttrium
    Zr = auto() # Zirconium
    Nb = auto() # Niobium
    Mo = auto() # Molybdenum
    Tc = auto() # Technetium
    Ru = auto() # Ruthenium
    Rh = auto() # Rhodium
    Pd = auto() # Palladium
    Ag = auto() # Silver
    Cd = auto() # Cadmium
    In = auto() # Indium
    Sn = auto() # Tin
    Sb = auto() # Antimony
    Te = auto() # Tellurium
    I  = auto() # Iodine
    Xe = auto() # Xenon
    Cs = auto() # Cesium
    Ba = auto() # Barium
    La = auto() # Lanthanum
    Ce = auto() # Cerium
    Pr = auto() # Praseodymium
    Nd = auto() # Neodymium
    Pm = auto() # Promethium
    Sm = auto() # Samarium    
    Eu = auto() # Europium
    Gd = auto() # Gadolinium
    Tb = auto() # Terbium
    Dy = auto() # Dysprosium
    Ho = auto() # Holmium
    Er = auto() # Erbium
    Tm = auto() # Thulium
    Yb = auto() # Ytterbium
    Lu = auto() # Lutetium
    Hf = auto() # Hafnium
    Ta = auto() # Tantalum
    W  = auto() # Tungsten
    Re = auto() # Rhenium
    Os = auto() # Osmium
    Ir = auto() # Iridium    
    Pt = auto() # Platinum
    Au = auto() # Gold
    Hg = auto() # Mercury 
    Tl = auto() # Thallium
    Pb = auto() # Lead 
    Bi = auto() # Bismuth
    Po = auto() # Polonium 
    At = auto() # Astatine
    Rn = auto() # Radon
    Fr = auto() # Francium
    Ra = auto() # Radium
    Ac = auto() # Actinium 
    Th = auto() # Thorium
    Pa = auto() # Protactinium
    U  = auto() # Uranium
    Np = auto() # Neptunium
    Pu = auto() # Plutonium
    Am = auto() # Americium
    Cm = auto() # Curium
    Bk = auto() # Berkelium
    Cf = auto() # Californium
    Es = auto() # Einsteinium
    Fm = auto() # Fermium
    Md = auto() # Mendelevium
    No = auto() # Nobelium
    Lr = auto() # Lawrencium
    Rf = auto() # Rutherfordium
    Db = auto() # Dubnium
    Sg = auto() # Seaborgium
    Bh = auto() # Bohrium
    Hs = auto() # Hassium
    Mt = auto() # Meitnerium
    Ds = auto() # Darmstadtium
    Rg = auto() # Roentgenium
    Cn = auto() # Copernicium
    Nh = auto() # Nihonium
    Fl = auto() # Flerovium
    Mc = auto() # Moscovium
    Lv = auto() # Livermorium
    Ts = auto() # Tennessine
    Og = auto() # Oganesson


DEFAULT_ATOM_RADII = {
    Symbol.H : 0.37,
    Symbol.He: 0.32,
    Symbol.Li: 1.28,
    Symbol.Be: 0.90,
    Symbol.B : 0.82,
    Symbol.C : 0.77,
    Symbol.N : 0.75,
    Symbol.O : 0.73,
    Symbol.F : 0.71,
    Symbol.Ne: 0.69,
    Symbol.Na: 1.54,
    Symbol.Mg: 1.30,
    Symbol.Al: 1.18,
    Symbol.Si: 1.11,
    Symbol.P : 1.06,
    Symbol.S : 1.02,
    Symbol.Cl: 0.99,
    Symbol.Ar: 0.97,
    Symbol.K : 2.03,
    Symbol.Ca: 1.74,
    Symbol.Sc: 1.44,
    Symbol.Ti: 1.32,
    Symbol.V : 1.22,
    Symbol.Cr: 1.18,
    Symbol.Mn: 1.17,
    Symbol.Fe: 1.17,
    Symbol.Co: 1.16,
    Symbol.Ni: 1.15,
    Symbol.Cu: 1.17,
    Symbol.Zn: 1.25,
    Symbol.Ga: 1.26,
    Symbol.Ge: 1.22,
    Symbol.As: 1.19,
    Symbol.Se: 1.20,
    Symbol.Br: 1.20,
    Symbol.Kr: 1.16,
    Symbol.Rb: 2.16,
    Symbol.Sr: 1.95,
    Symbol.Y : 1.88,
    Symbol.Zr: 1.75,
    Symbol.Nb: 1.64,
    Symbol.Mo: 1.54,
    Symbol.Tc: 1.47,
    Symbol.Ru: 1.46,
    Symbol.Rh: 1.42,
    Symbol.Pd: 1.39,
    Symbol.Ag: 1.45,
    Symbol.Cd: 1.48,
    Symbol.In: 1.42,
    Symbol.Sn: 1.39,
    Symbol.Sb: 1.39,
    Symbol.Te: 1.38,
    Symbol.I : 1.39,
    Symbol.Xe: 1.40,
    Symbol.Cs: 2.44,
    Symbol.Ba: 2.15,
    Symbol.La: 2.07,
    Symbol.Ce: 2.04,
    Symbol.Pr: 2.03,
    Symbol.Nd: 2.01,
    Symbol.Pm: 1.99,
    Symbol.Sm: 1.98,
    Symbol.Eu: 1.98,
    Symbol.Gd: 1.96,
    Symbol.Tb: 1.94,
    Symbol.Dy: 1.92,
    Symbol.Ho: 1.92,
    Symbol.Er: 1.89,
    Symbol.Tm: 1.90,
    Symbol.Yb: 1.87,
    Symbol.Lu: 1.87,
    Symbol.Hf: 1.75,
    Symbol.Ta: 1.70,
    Symbol.W : 1.62,
    Symbol.Re: 1.51,
    Symbol.Os: 1.44,
    Symbol.Ir: 1.41,
    Symbol.Pt: 1.36,
    Symbol.Au: 1.36,
    Symbol.Hg: 1.32,
    Symbol.Tl: 1.45,
    Symbol.Pb: 1.46,
    Symbol.Bi: 1.48,
    Symbol.Po: 1.40,
    Symbol.At: 1.50,
    Symbol.Rn: 1.50,
    Symbol.Fr: 2.60,
    Symbol.Ra: 2.21,
    Symbol.Ac: 2.15,
    Symbol.Th: 2.06,
    Symbol.Pa: 2.00,
    Symbol.U : 1.96,
    Symbol.Np: 1.90,
    Symbol.Pu: 1.87,
    Symbol.Am: 1.80,
    Symbol.Cm: 1.69,
    Symbol.Bk: 1.66,
    Symbol.Cf: 1.65,
    Symbol.Es: 1.60,
    Symbol.Fm: 1.57,
    Symbol.Md: 1.56,
    Symbol.No: 1.54,
    Symbol.Lr: 1.56,
    Symbol.Rf: 1.57,
    Symbol.Db: 1.57,
    Symbol.Sg: 1.57,
    Symbol.Bh: 1.57,
    Symbol.Hs: 1.57,
    Symbol.Mt: 1.57,
    Symbol.Ds: 1.57,
    Symbol.Rg: 1.57,
    Symbol.Cn: 1.57,
    Symbol.Nh: 1.57,
    Symbol.Fl: 1.57,
    Symbol.Mc: 1.57,
    Symbol.Lv: 1.57,
    Symbol.Ts: 1.57,
    Symbol.Og: 1.57
}


def default_atom_radius(symbol: Symbol, factor: float = 1.0) -> float:
    """
    Returns the default atom radius for the given symbol.

    Arguments
    ---------
    symbol : Symbol
        The symbol of the atom.
    factor : float
        The factor to scale the radius.
    
    Returns
    -------
    float
        The default atom radius.
    """
    return DEFAULT_ATOM_RADII[symbol] * factor


def default_atom_color(symbol: Symbol) -> Color:
    """
    Returns the default atom color for the given symbol.

    Arguments
    ---------
    symbol : Symbol
        The symbol of the atom.
    
    Returns
    -------
    Color
        The default atom color.
    """
    match symbol:
        case Symbol.H : color = Color(255/255, 255/255, 255/255)
        case Symbol.C : color = Color(211/255, 211/255, 211/255)
        case Symbol.N : color = Color(  0/255,   0/255, 255/255)
        case Symbol.O : color = Color(255/255,   0/255,   0/255)
        case Symbol.S : color = Color(255/255, 255/255,   0/255)
        case Symbol.P : color = Color(255/255, 128/255,   0/255)
        case Symbol.F : color = Color(  0/255, 255/255,   0/255)
        case Symbol.Cl: color = Color(  0/255, 255/255,   0/255)
        case Symbol.Br: color = Color(  0/255, 255/255,   0/255)
        case Symbol.I : color = Color(  0/255, 255/255, 255/255)
        case _        : color = Color(  0/255,   0/255,   0/255)

    return color


@dataclass
class Atom:
    """
    An atom describes a node in a molecular graph.

    Arguments
    ---------
    symbol : Symbol
        The symbol of the atom.
    position : Point
        The position of the atom.
    radius : float
        The radius of the atom.
    color : Color
        The color of the atom.
    theta_resolution : int
        The number of points used to draw the atom.
    phi_resolution : int
        The number of points used to draw the atom.
    """
    symbol: Symbol
    position: Point

    radius: Radius = default_atom_radius(Symbol)
    color: Color = default_atom_color(Symbol)

    theta_resolution: int = 90
    phi_resolution: int = 90


# =============================================================================
# Bond definitions and functionalities.
# =============================================================================
DEFAULT_BOND_RADIUS = 0.1


class BondType(Enum):
    """
    Defines the type of a bond.
    """
    Single = auto()
    Double = auto()
    Triple = auto()
    Aromatic = auto()


@dataclass 
class Bond:
    """
    A bond describes an edge in a molecular graph.

    Arguments
    ---------
    type : BondType
        The type of the bond.
    segment_a : Point
        Position one end of bond.
    segment_b : Point
        Position other end of bond.
    color_segment_a : Color
        Color of segment one end of bond.
    color_segment_b : Color
        Color of segment other end of bond.
    radius : float
        The radius of the bond.
    resolution : int 
        The number of points used to draw the bond.
    n_sides : int 
        The number of sides used to draw the bond.
    """
    type: BondType

    segment_a: Point 
    segment_b: Point 
    color_segment_a: Color 
    color_segment_b: Color
    radius: Radius = DEFAULT_BOND_RADIUS

    resolution: int = 100
    n_sides: int = 100


# =============================================================================
# Molecule definitions and functionalities.
# =============================================================================
@dataclass
class Molecule:
    """
    A molecule is a collection of atoms and bonds.
    """
    atoms: List[Atom]
    bonds: List[Bond]


# =============================================================================
# Geometry utilities.
# =============================================================================
def centroid(arr: np.array) -> np.array:
    """
    Calculates centroid of x, y, z coordinates of point cloud. 

    Arguments
    ---------
    arr : np.array
        Array of x, y, z coordinates with dimensions (N, 3).

    Returns
    -------
    centroid : np.array
        Centroid of point cloud with dimensions (1, 3). 
    """
    return sum(arr, axis=0) / len(arr)


def rotation_matrix_x(r: float) -> np.array:
    """
    Calculates rotation matrix for rotation around x-axis.

    Arguments
    ---------
    r : float
        Rotation angle in radians.

    Returns
    -------
    rot_matrix : np.array
        Rotation matrix for rotation around x-axis.

    For more information on 4x4 rotation matrices:
    https://iwatobipen.wordpress.com/author/iwatobipen/ 
    http://www.cs.cmu.edu/afs/cs/academic/class/15462-s10/www/lec-slides/lec04.pdf
    """
    return np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, np.cos(r), -np.sin(r), 0.0],
        [0.0, np.sin(r), np.cos(r), 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])


def rotation_matrix_y(r: float) -> np.array:
    """
    Calculates rotation matrix for rotation around y-axis.
    
    Arguments
    ---------
    r : float 
        Rotation angle in radians.
    
    Returns
    -------
    rot_matrix : np.array
        Rotation matrix for rotation around y-axis.

    For more information on 4x4 rotation matrices:
    https://iwatobipen.wordpress.com/author/iwatobipen/ 
    http://www.cs.cmu.edu/afs/cs/academic/class/15462-s10/www/lec-slides/lec04.pdf
    """
    return np.array([
        [np.cos(r), 0.0, np.sin(r), 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [-np.sin(r), 0.0, np.cos(r), 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])


def rotation_matrix_z(r: float) -> np.array:
    """
    Calculates rotation matrix for rotation around z-axis.
    
    Arguments
    ---------
    r : float
        Rotation angle in radians.
    
    Returns
    -------
    rot_matrix : np.array
        Rotation matrix for rotation around z-axis.

    For more information on 4x4 rotation matrices:
    https://iwatobipen.wordpress.com/author/iwatobipen/ 
    http://www.cs.cmu.edu/afs/cs/academic/class/15462-s10/www/lec-slides/lec04.pdf
    """
    return np.array([
        [np.cos(r), -np.sin(r), 0.0, 0.0],
        [np.sin(r), np.cos(r), 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])


# =============================================================================
# Main functionality for drawing.
# =============================================================================
def create_atom_mesh(atom: Atom) -> Sphere:
    """
    Creates a mesh for an atom.

    Arguments
    ---------
    atom : Atom
        The atom to create the mesh for.

    Returns
    -------
    atom_mesh : Sphere
        The mesh for the atom.
    """
    return Sphere(
        center=atom.position.as_array,
        direction=atom.position.as_array,
        radius=atom.radius,
        theta_resolution=atom.theta_resolution,
        phi_resolution=atom.phi_resolution,
    )


def create_bond_mesh(bond: Bond) -> Tuple[Tube, Tube]:
    """
    Creates a mesh for a bond.

    Arguments
    ---------
    bond : Bond
        The bond to create the mesh for.

    Returns
    -------
    Tuple[Tube, Tube]
        The meshes for both bond segments.
    """
    bond_center = centroid([bond.segment_a.as_array, bond.segment_b.as_array])

    segment_a = Tube(
        pointa=bond.segment_a.as_array,
        pointb=bond_center,
        resolution=bond.resolution,
        n_sides=bond.n_sides,
        radius=bond.radius,
    )
    segment_b = Tube(
        pointa=bond.segment_b.as_array,
        pointb=bond_center,
        resolution=bond.resolution,
        n_sides=bond.n_sides,
        radius=bond.radius,
    )

    return segment_a, segment_b


@dataclass
class Drawer:
    """
    Defines the drawer for the scene.

    Arguments
    ---------
    background_color : Color
        Background color of the scene.
    off_screen : bool
        If True, the scene is drawn on an off-screen buffer.
    """
    background_color: Color = Color(255/255, 255/255, 255/255)
    off_screen: bool = True 

    plotter = Plotter(off_screen=off_screen)
    plotter.background_color = background_color

    def add_molecule(self, molecule: Molecule) -> None:
        """
        Adds a molecule to the scene.

        Arguments
        ---------
        molecule : Molecule
            Molecule to be added to the scene.
        """
        # Draw atoms in scene.
        for atom in molecule.atoms: 
            atom_mesh = create_atom_mesh(atom)
            self.plotter.add_mesh(atom_mesh, color=atom.color)
        
        # Draw bonds in scene.
        for bond in molecule.bonds:
            bond_mesh_segment_a, bond_mesh_segment_b = create_bond_mesh(bond)
            self.plotter.add_mesh(bond_mesh_segment_a, color=bond.color_segment_a)
            self.plotter.add_mesh(bond_mesh_segment_b, color=bond.color_segment_b)

    def show(self) -> None:
        """
        Draws the scene.
        """
        self.plotter.show()

    @error_handler
    def save_fig(
        self,
        filename: str,
        transparent_background: bool = False,
        window_size: Tuple[int, int] = (800, 800)
    ) -> None:
        """
        Draws the scene and saves it to output file.

        Arguments
        ---------
        filename : str
            Name of the output file.
        transparent_background : bool
            If True, the background is transparent.
        window_size : Tuple[int, int]
            Size of the window.
        """
        self.plotter.screenshot(
            filename, 
            transparent_background=transparent_background,
            window_size=window_size
        )
        