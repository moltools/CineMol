[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Maintainer](https://img.shields.io/badge/Maintainer-davidmeijer-blue)](https://github.com/davidmeijer)
[![Generic badge](https://img.shields.io/badge/Version-alpha-green.svg)](https://shields.io/)

# CineMol

<img src="./logo.png" alt="logo" width="100">

CineMol is a direct-to-SVG small molecule drawer. 

Try it out with a low resolution [here](https://moltools.nl/cinemol).

## Usage

CineMol will be made available as command line tool and a NuGet package for .NET 6.0 shortly.

See `scripts/draw_molecule.fsx` for an example script which you can use right now.

You can execute this script yourself. First build CineMol by running on the command line:

```bash
dotnet build src/CineMol/CineMol.fsproj
```

Then run the script with:

```bash 
dotnet fsi scripts/draw_molecule.fsx <path_to_input_sdf_file> <path_to_output_svg_file>
```

or 

```bash
dotnet fsi scripts/draw_molecule.fsx data/penicillin_G.sdf output/penicillin_G.svg
```

The SDF file for penicillin G referenced above was downloaded from [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/Penicillin-G).

The output SVG file will look like this:

<table>
  <tr>
    <td>No rotation</td>
    <td>90&deg; rotation over X-axis</td>
    <td>beta-lactam ring highlighted</td>
  </tr>
  <tr>
    <td>
        <img src="./data/penicillin_G.svg" alt="example1" width="160">
    </td>
    <td>
        <img src="./data/penicillin_G_rotated.svg" alt="example2" width="160">
    </td>
    <td>
        <img src="./data/penicillin_G_highlighted.svg" alt="example3" width="160">
    </td>
  </tr>
 </table>

## Interactive viewer 

For an explanation on how to start an interactive molecule viewer locally, please see the `app` folder.

## Using CineMol with Python

You can easily use [Fable Python](https://fable.io/Fable.Python/) in order to use CineMol from Python.


