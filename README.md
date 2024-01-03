[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# CineMol

<img src="./logo.png" alt="logo" width="100">

CineMol is a direct-to-SVG small molecule drawer. 

You can install CineMol with pip from the root of this repository:

```bash
pip install .
```

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
 </table>

## Usage Python

```python
from cinemol.parsers import parse_sdf 
from cinemol.chemistry import Style, Look, draw_molecule

with open("path/to/your/molecule.sdf", "r") as f:
    sdf_str = f.read()

atoms, bonds = parse_sdf(sdf_str)
svg_str = draw_molecule(atoms, bonds, style=Style.Tube, look=Look.Glossy, resolution=100, scale=10.0)
```

See `src/cinemol/chemistry.py` for more options.

## Usage CLI

```bash
cinemol -i path/to/your/molecule.sdf -o path/to/your/molecule.svg
```

Run `cinemol -h` for more options.
