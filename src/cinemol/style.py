# -*- coding: utf-8 -*-

"""This module contains classes for styling molecular depictions."""

import typing as ty
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from cinemol.geometry import Point2D

# ==============================================================================
# Color
# ==============================================================================


class Color:
    """A color in RGB format."""

    def __init__(self, r: int, g: int, b: int) -> None:
        """Initialize the color.

        :param r: The red component of the color.
        :type r: int
        :param g: The green component of the color.
        :type g: int
        :param b: The blue component of the color.
        :type b: int
        :raises TypeError: If r, g, or b is not an integer.
        """
        if not isinstance(r, int) or not isinstance(g, int) or not isinstance(b, int):
            raise TypeError(
                f"Expected (r=int, g=int, b=int), got'(r={type(r)}, g={type(g)}, b={type(b)})'"
            )

        self.r = max(0, min(255, r))
        self.g = max(0, min(255, g))
        self.b = max(0, min(255, b))

    def __str__(self) -> str:
        """Return the string representation of the color.

        :return: The string representation of the color.
        :rtype: str
        """
        return f"rgb({self.r},{self.g},{self.b})"

    def __eq__(self, other: ty.Any) -> bool:
        """Check if self is equal to other.

        :param other: Another object.
        :type other: ty.Any
        :return: If self is equal to other.
        :rtype: bool
        """
        if not isinstance(other, Color):
            return False

        return (self.r == other.r) and (self.g == other.g) and (self.b == other.b)

    def to_hex(self) -> str:
        """Convert color to a hex string.

        :return: Hex representation of the color.
        :rtype: str
        """
        return f"#{self.r:02x}{self.g:02x}{self.b:02x}"

    def diffuse(self, alpha: float) -> "Color":
        """Diffuse color with alpha value.

        :param alpha: Diffusion factor.
        :type alpha: float
        :return: Diffused color.
        :rtype: Color
        :raises TypeError: If alpha is not a float.
        """
        if not isinstance(alpha, float):
            raise TypeError(f"Expected alpha to be float, got '{type(alpha)}'")

        alpha = max(0, min(1, alpha))
        color = Color(int(self.r * alpha), int(self.g * alpha), int(self.b * alpha))

        return color


# ==============================================================================
# Atom colors
# ==============================================================================


class AtomColoringScheme(ABC):
    """Abstract base class for atom coloring schemes."""

    @abstractmethod
    def get_color(self, atom_symbol: str) -> Color:
        """
        Return the color of an atom.

        :param atom_symbol: The symbol of the atom.
        :type atom_symbol: str
        :return: The color of the atom.
        :rtype: Color
        """


class CoreyPaulingKoltungAtomColor(AtomColoringScheme):
    """Corey-Pauling-Koltun coloring convention for atoms.

    Source: https://en.wikipedia.org/wiki/CPK_coloring
    """

    H = Color(255, 255, 255)  # White
    C = Color(80, 80, 80)  # Dark gray
    N = Color(0, 0, 255)  # Blue
    Oxygen = Color(255, 0, 0)  # Red
    P = Color(255, 165, 0)  # Orange
    S = Color(255, 255, 0)  # Yellow
    B = Color(245, 245, 220)  # Beige
    Br = Color(139, 0, 0)  # Dark red
    Iodine = Color(148, 0, 211)  # Dark violet
    Ti = Color(128, 128, 128)  # Gray
    Fe = Color(255, 140, 0)  # Dark orange
    F = Color(0, 128, 0)  # Green
    Cl = Color(0, 128, 0)  # Green
    He = Color(0, 255, 255)  # Cyan
    Ne = Color(0, 255, 255)  # Cyan
    Ar = Color(0, 255, 255)  # Cyan
    Kr = Color(0, 255, 255)  # Cyan
    Xe = Color(0, 255, 255)  # Cyan
    Li = Color(238, 130, 238)  # Violet
    Na = Color(238, 130, 238)  # Violet
    K = Color(238, 130, 238)  # Violet
    Rb = Color(238, 130, 238)  # Violet
    Cs = Color(238, 130, 238)  # Violet
    Fr = Color(238, 130, 238)  # Violet
    Be = Color(0, 100, 0)  # Dark green
    Mg = Color(0, 100, 0)  # Dark green
    Ca = Color(0, 100, 0)  # Dark green
    Sr = Color(0, 100, 0)  # Dark green
    Ba = Color(0, 100, 0)  # Dark green
    Ra = Color(0, 100, 0)  # Dark green
    Cd = Color(170, 51, 106)  # Dark pink

    def get_color(self, atom_symbol: str) -> Color:
        """
        Return the color of an atom.

        :param atom_symbol: The symbol of the atom.
        :type atom_symbol: str
        :return: The color of the atom.
        :rtype: Color
        """
        default = Color(255, 192, 203)  # Pink

        # We define O and I separately because otherwise the variable name would
        # be ambiguous according to PEP8.
        if atom_symbol == "O":
            return self.Oxygen
        elif atom_symbol == "I":
            return self.Iodine
        else:
            return getattr(self, atom_symbol, default)


# ==============================================================================
# Atom radius
# ==============================================================================


class AtomRadiusScheme(ABC):
    """Abstract base class for atom radius schemes."""

    @abstractmethod
    def to_angstrom(self, atom_symbol: str) -> float:
        """Return the radius of an atom in angstrom.

        :param atom_symbol: The symbol of the atom.
        :type atom_symbol: str
        :return: The radius of the atom in angstrom.
        :rtype: float
        """


class PubChemAtomRadius(AtomRadiusScheme):
    """Atomic radii (van der Waals) in picometer from PubChem.

    Source: https://pubchem.ncbi.nlm.nih.gov/periodic-table/#property=AtomicRadius
    """

    H = 120.0
    He = 140.0
    Li = 182.0
    Be = 153.0
    B = 192.0
    C = 170.0
    N = 155.0
    Oxygen = 152.0
    F = 135.0
    Ne = 154.0
    Na = 227.0
    Mg = 173.0
    Al = 184.0
    Si = 210.0
    P = 180.0
    S = 180.0
    Cl = 175.0
    Ar = 188.0
    K = 275.0
    Ca = 231.0
    Sc = 211.0
    Ti = 187.0
    V = 179.0
    Cr = 189.0
    Mn = 197.0
    Fe = 194.0
    Co = 192.0
    Ni = 163.0
    Cu = 140.0
    Zn = 139.0
    Ga = 187.0
    Ge = 211.0
    As = 185.0
    Se = 190.0
    Br = 183.0
    Kr = 202.0
    Rb = 303.0
    Sr = 249.0
    Y = 219.0
    Zr = 186.0
    Nb = 207.0
    Mo = 209.0
    Tc = 209.0
    Ru = 207.0
    Rh = 195.0
    Pd = 202.0
    Ag = 172.0
    Cd = 158.0
    In = 193.0
    Sn = 217.0
    Sb = 206.0
    Te = 206.0
    Iodine = 198.0
    Xe = 216.0
    Cs = 343.0
    Ba = 268.0
    Lu = 221.0
    Hf = 212.0
    Ta = 217.0
    W = 210.0
    Re = 217.0
    Os = 216.0
    Ir = 202.0
    Pt = 209.0
    Au = 166.0
    Hg = 209.0
    Tl = 196.0
    Pb = 202.0
    Bi = 207.0
    Po = 197.0
    At = 202.0
    Rn = 220.0
    Fr = 348.0
    Ra = 283.0

    def to_angstrom(self, atom_symbol: str) -> float:
        """Return the radius of an atom in angstrom.

        :param atom_symbol: The symbol of the atom.
        :type atom_symbol: str
        :return: The radius of the atom in angstrom.
        :rtype: float
        """
        default = 170.0
        factor = 0.01  # Convert from picometer to angstrom.

        # We define O and I separately because otherwise the variable name would
        # be ambiguous according to PEP8.
        if atom_symbol == "O":
            return self.Oxygen * factor
        elif atom_symbol == "I":
            return self.Iodine * factor
        else:
            return getattr(self, atom_symbol, default) * factor


# ==============================================================================
# Art styles
# ==============================================================================


class Depiction:
    """Abstract base class for depiction styles."""


@dataclass
class Cartoon(Depiction):
    """Cartoon fill style.

    :param fill_color: The color of the fill.
    :type fill_color: Color
    :param outline_color: The color of the stroke.
    :type outline_color: Color
    :param outline_width: The width of the stroke.
    :type outline_width: float
    :param opacity: The opacity of the fill.
    :type opacity: float
    """

    fill_color: Color
    outline_color: Color = field(default_factory=lambda: Color(0, 0, 0))
    outline_width: float = 0.05
    opacity: float = 1.0


@dataclass
class Glossy(Depiction):
    """Glossy fill style.

    :param fill_color: The color of the fill.
    :param outline_color: The color of the stroke.
    """

    fill_color: Color
    opacity: float = 1.0


class FillStyle:
    """Abstract base class for fill styles."""


class Wire(FillStyle):
    """Wire fill style."""

    def __init__(self, stroke_color: Color, stroke_width: float, opacity: float) -> None:
        """Initialize the wire fill style.

        :param stroke_color: The color of the stroke.
        :type stroke_color: Color
        :param stroke_width: The width of the stroke.
        :type stroke_width: float
        :param opacity: The opacity of the stroke.
        :type opacity: float
        """
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
        self.opacity = opacity


class Solid(FillStyle):
    """Solid fill style."""

    def __init__(
        self, fill_color: Color, stroke_color: Color, stroke_width: float, opacity: float
    ) -> None:
        """Initialize the solid fill style.

        :param fill_color: The color of the fill.
        :type fill_color: Color
        :param stroke_color: The color of the stroke.
        :type stroke_color: Color
        :param stroke_width: The width of the stroke.
        :type stroke_width: float
        :param opacity: The opacity of the fill.
        :type opacity: float
        """
        self.fill_color = fill_color
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
        self.opacity = opacity


class RadialGradient(FillStyle):
    """Radial gradient fill style."""

    def __init__(self, fill_color: Color, center: Point2D, radius: float, opacity: float) -> None:
        """Initialize the radial gradient fill style.

        :param fill_color: The color of the radial gradient.
        :type fill_color: Color
        :param center: The center of the radial gradient.
        :type center: Point2D
        :param radius: The radius of the radial gradient.
        :type radius: float
        :param opacity: The opacity of the radial gradient.
        :type opacity: float
        """
        self.fill_color = fill_color
        self.center = center
        self.radius = radius
        self.opacity = opacity


class LinearGradient(FillStyle):
    """Linear gradient fill style."""

    def __init__(self, fill_color: Color, start: Point2D, end: Point2D, opacity: float) -> None:
        """Initialize the linear gradient fill style.

        :param fill_color: The color of the linear gradient.
        :type fill_color: Color
        :param start: The start of the linear gradient.
        :type start: Point2D
        :param end: The end of the linear gradient.
        :type end: Point2D
        :param opacity: The opacity of the linear gradient.
        :type opacity: float
        """
        self.fill_color = fill_color
        self.start = start
        self.end = end
        self.opacity = opacity


@dataclass
class Fill:
    """Fill style.

    :param reference: The reference of the fill style.
    :type reference: str
    :param fill_style: The fill style.
    :type fill_style: FillStyle
    """

    reference: str
    fill_style: FillStyle

    def to_svg(self) -> ty.Tuple[str, ty.Optional[str]]:
        """Return the SVG representation of the fill style.

        :return: The SVG representation of the fill style. The first string is the
            style string, the second string is the definition string.
        :rtype: str, str
        :raises TypeError: If the fill style is not a FillStyle.
        """
        if isinstance(self.fill_style, Wire):
            stroke_color = self.fill_style.stroke_color.to_hex()
            stroke_width = self.fill_style.stroke_width
            opacity = self.fill_style.opacity

            style_str = (
                f".{self.reference}"
                f"{{stroke:{stroke_color};"
                f"stroke-width:{stroke_width:.3f}px;"
                f"stroke-opacity:{opacity};"
                f"stroke-linecap:round;"
                f"stroke-linejoin:round;}}"
            )
            definition_str = None

            return style_str, definition_str

        elif isinstance(self.fill_style, Solid):
            fill_color = self.fill_style.fill_color.to_hex()
            stroke_color = self.fill_style.stroke_color.to_hex()
            stroke_width = self.fill_style.stroke_width
            opacity = self.fill_style.opacity

            style_str = (
                f".{self.reference}"
                f"{{fill:{fill_color};"
                f"stroke:{stroke_color};"
                f"stroke-width:{stroke_width:.3f}px;"
                f"opacity:{opacity};}}"
            )
            definition_str = None

            return style_str, definition_str

        elif isinstance(self.fill_style, RadialGradient):
            cx = self.fill_style.center.x
            cy = self.fill_style.center.y
            r = self.fill_style.radius

            stop_color_offset_a = self.fill_style.fill_color.to_hex()
            stop_color_offset_b = self.fill_style.fill_color.diffuse(0.75).to_hex()

            style_str = f".{self.reference}{{fill:url(#{self.reference});}}"
            definition_str = (
                "<radialGradient"
                f' id="{self.reference}"'
                f' cx="{cx:.3f}" cy="{cy:.3f}"'
                f' r="{r:.3f}" fx="{cx:.3f}" fy="{cy:.3f}"'
                ' gradientTransform="matrix(1,0,0,1,0,0)"'
                ' gradientUnits="userSpaceOnUse"'
                f' opacity="{self.fill_style.opacity}">'
                f'<stop offset="0.00" stop-color="{stop_color_offset_a}"/>'
                f'<stop offset="1.00" stop-color="{stop_color_offset_b}"/>'
                "</radialGradient>"
            )

            return style_str, definition_str

        elif isinstance(self.fill_style, LinearGradient):
            multiplier = 3.0  # Makes linear gradient look better.
            x1 = self.fill_style.start.x * multiplier
            y1 = self.fill_style.start.y * multiplier
            x2 = self.fill_style.end.x * multiplier
            y2 = self.fill_style.end.y * multiplier

            stop_color_offset_a = self.fill_style.fill_color.to_hex()
            stop_color_offset_b = self.fill_style.fill_color.diffuse(0.75).to_hex()

            style_str = f".{self.reference}{{fill:url(#{self.reference});}}"
            definition_str = (
                "<linearGradient"
                f' id="{self.reference}"'
                f' x1="{x1:.3f}" y1="{y1:.3f}"'
                f' x2="{x2:.3f}" y2="{y2:.3f}"'
                ' gradientUnits="userSpaceOnUse"'
                ' spreadMethod="reflect"'
                ' gradientTransform="rotate(90)"'
                f' opacity="{self.fill_style.opacity}">'
                f'<stop offset="0.00" stop-color="{stop_color_offset_a}"/>'
                f'<stop offset="1.00" stop-color="{stop_color_offset_b}"/>'
                "</linearGradient>"
            )

            return style_str, definition_str

        else:
            raise TypeError(f"Expected FillStyle, got '{type(self.fill_style)}'")
