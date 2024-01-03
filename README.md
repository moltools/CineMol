[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# CineMol

<img src="./logo.png" alt="logo" width="100">

CineMol is a direct-to-SVG small molecule drawer. 

## Installation

You can install CineMol with pip from the root of this repository:

```bash
pip install .
```

## Usage CLI

```bash
cinemol -i path/to/your/molecule.sdf -o path/to/your/molecule.svg -s tube -l glossy -r 100 -sc 10.0 --hs
```

Command line options:

* `-i`: input file path to SDF file (str).
* `-o`: output file path to SVG file (str).
* `-s`: style (`spacefilling`/`ballandstick`/`tube`/`wireframe`) (str).
* `-l`: look (`cartoon`/`glossy`) (str).
* `-r`: resolution of the SVG that determines the number of points to be drawn on the sphere and cylinder geometries in order to approximate them (int).
* `-sc`: scales radii, coordinates, and strokes by this factor (float).
* `-rx`: rotation over the x-axis in degrees (float).
* `-ry`: rotation over the y-axis in degrees (float).
* `-rz`: rotation over the z-axis in degrees (float).
* `--hs`: show hydrogens (bool).
* `--vb`: verbose (bool).
* `-v`: print version to stdout.

### Styling options

A penicillin G conformer, [retrieved from PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/Penicillin-G), drawn with CineMol in various looks and styles:

<table>
  <tr>
    <td>Space-filling cartoon</td>
    <td>Space-filling glossy</td>
  </tr>
  <tr>
    <td><img src="svgs/cartoon_spacefilling.svg" width=200 height=125></td>
    <td><img src="svgs/glossy_spacefilling.svg" width=200 height=125></td>
  </tr>
  <tr>
    <td>Ball-and-stick cartoon</td>
    <td>Ball-and-stick glossy</td>
  </tr>
  <tr>
    <td><img src="svgs/cartoon_ballandstick.svg" width=200 height=125></td>
    <td><img src="svgs/glossy_ballandstick.svg" width=200 height=125></td>
  </tr>
  <tr>
    <td>Tube cartoon</td>
    <td>Tube glossy</td>
  </tr>
  <tr>
    <td><img src="svgs/cartoon_tube.svg" width=200 height=125></td>
    <td><img src="svgs/glossy_tube.svg" width=200 height=125></td>
  </tr>
  <tr>
    <td colspan="2">Wireframe</td>
  </tr>
  <tr>
    <td colspan="2"><img src="svgs/wireframe.svg" width=400 height=125></td>
  </tr>
 </table>

## Usage Python

```python
from cinemol.parsers import parse_sdf 
from cinemol.chemistry import Style, Look, draw_molecule

# Parse atoms and bonds from an SDF file ...
with open("path/to/your/molecule.sdf", "r") as f:
    sdf_str = f.read()

atoms, bonds = parse_sdf(sdf_str)

# ... or create your own atom and bond objects:
atoms = [
    Atom(0, "C", (0.0, 0.0, 0.0)), 
    Atom(1, "N", (1.0, 0.0, 0.0))
]
bonds = [
    Bond(0, 1, order=3)
]

# Draw molecule to SVG string:
svg_str = draw_molecule(atoms, bonds, style=Style.Tube, look=Look.Glossy, resolution=100, scale=10.0)
```

Options for `draw_molecule`:

* `atoms`: list of `Atom` objects
* `bonds`: list of `Bond` objects
* `style`: `Style` enum (`SpaceFilling`/`BallAndStick`/`Tube`/`Wireframe`)
* `look`: `Look` enum (`Cartoon`/`Glossy`)
* `resolution`: resolution of the SVG that determines the number of points to be drawn on the sphere and cylinder geometries in order to approximate them (int).
* `rotation_over_x_axis`: rotation over the x-axis in degrees (float).
* `rotation_over_y_axis`: rotation over the y-axis in degrees (float).
* `rotation_over_z_axis`: rotation over the z-axis in degrees (float).
* `verbose`: print drawing progress to stdout (bool).
* `scale`: scales radii, coordinates, and strokes by this factor (float).

Initialize an `Atom` object:

* `index`: index of the atom in the molecule (int). Used to look up bond coordinates.
* `symbol`: element symbol (str).
* `coordinates`: 3D coordinates of the atom (enumerable of three floats as x, y, z coordinates).
* `radius`: radius of the atom (float). If not supplied, PubChem radius is used based on atom symbol (for PuBchem radii, see `src/cinemol/style.py`).
* `color`: color of atom (enumerable of three ints as r, g, b values). If not supplied, CPK color is used based on atom symbol (for CPK colors, see `src/cinemol/style.py`).
* `opacity`: opacity of atom (float between 0.0 and 1.0). If not supplied, 1.0 is used.

Initialize a `Bond` object:

* `start_index`: index of the start atom in the molecule (int). Used to look up bond coordinates.
* `end_index`: index of the end atom in the molecule (int). Used to look up bond coordinates.
* `order`: bond order (int). If not supplied, 1 is used.
* `color`: color of bond (enumerable of three ints as r, g, b values). If not supplied, color is used from nearest atom.
* `radius`: radius of bond (float). If not supplied, 0.2 is used. 
* `opacity`: opacity of bond (float between 0.0 and 1.0). If not supplied, 1.0 is used.

### Fine-grained control

See `Scene` in `src/cinemol/model.py` if you want direct access to the drawing scene. The scene can be giving various objects for drawing (e.g., `ModelSphere`, `ModelCylinder` ,`ModelWire`). The `draw_molecule` function is a convenience function on top of `Scene` that creates a scene, adds the molecule to the scene, and draws the scene to an SVG string. 