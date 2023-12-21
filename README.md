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

## Interactive viewer 

For an explanation on how to start an interactive molecule viewer locally, please see the `app` folder.

## Using CineMol with Python

You can easily use [Fable Python](https://fable.io/Fable.Python/) in order to use CineMol from Python (>=3.10):

```bash 
dotnet new tool-manifest
dotnet tool install fable --prerelease
dotnet tool run fable --lang Python CineMol.fsproj
```

Since CineMol-py has already been created (see: `/src/cinemol_py`) you can install it directly with:

```bash
pip install .
```

or by creating a dedicated environment first, together with the RDKit cheminformatics toolkit:

```bash
conda env create -f environment.yml
```




