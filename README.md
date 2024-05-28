<p align="center">
  <img src="https://github.com/moltools/CineMol/blob/main/logo.png" height="150">
</p>

<h1 align="center">
  CineMol
</h1>

<p align="center">
    <a href="https://github.com/MolTools/CineMol/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/MolTools/CineMol/actions/workflows/tests.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/cinemol">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/cinemol" /></a>
    <a href="https://pypi.org/project/cinemol">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/cinemol" /></a>
    <a href="https://github.com/MolTools/CineMol/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/cinemol" /></a>
    <a href='https://cinemol.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/cinemol/badge/?version=latest' alt='Documentation Status' /></a>
    <a href="https://codecov.io/gh/MolTools/CineMol/branch/main">
        <img src="https://codecov.io/gh/MolTools/CineMol/branch/main/graph/badge.svg" alt="Codecov status" /></a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /></a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' /></a>
    <a href="https://github.com/MolTools/CineMol/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/></a>
    <a href="https://doi.org/10.5281/zenodo.11242217">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.11242217.svg" alt="DOI"></a>
</p>

CineMol is a direct-to-SVG small molecule drawer. Read more about CineMol in our [publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00851-y). 

You can try out CineMol online [here](https://moltools.bioinformatics.nl/cinemol).


## 💪 Getting Started

The cinemol command line tool is automatically installed. It can
be used from the shell with the `--help` flag to show all subcommands:

```shell
python3 -m cinemol --help
```

The `cinemol` command line tool can be used to convert 3D molecular structures to SVG images. The following example
converts a penicillin G conformer to a spacefilling SVG image:

```shell
python3 -m cinemol examples/data/penicillin_G.sdf examples/svgs/cartoon_spacefilling.svg -s spacefilling -l cartoon -r 100 -sc 10.0 -hs
```

This command will generate the following SVG image:

<p align="center">
    <img src="https://github.com/moltools/CineMol/blob/main/examples/svgs/cartoon_spacefilling.svg" width=200 height=125>
</p>

CineMol can also be used as a Python package. You can find examples in the [examples](examples) folder.

Find more information in the [documentation](https://cinemol.readthedocs.io/en/latest/).

## 🚀 Installation

The most recent release can be installed from
[PyPI](https://pypi.org/project/cinemol/) with:

```shell
pip install cinemol
```

The most recent code and data can be installed directly from GitHub with:

```shell
pip install git+https://github.com/MolTools/CineMol.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/MolTools/CineMol/blob/main/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.

### 📖 Citation

If you use CineMol in your research, please [cite](CITATION.cff) our pre-print.

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
git clone git+https://github.com/MolTools/CineMol.git
cd CineMol
pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
tox
```

Additionally, these tests are automatically re-run with each commit in a
[GitHub Action](https://github.com/MolTools/CineMol/actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```shell
git clone git+https://github.com/MolTools/CineMol.git
cd CineMol
tox -e docs
open docs/build/html/index.html
``` 

The documentation automatically installs the package as well as the `docs`
extra specified in the [`setup.cfg`](setup.cfg). `sphinx` plugins
like `texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

The documentation can be deployed to [ReadTheDocs](https://readthedocs.io) using 
[this guide](https://docs.readthedocs.io/en/stable/intro/import-guide.html).
The [`.readthedocs.yml`](.readthedocs.yml) YAML file contains all the configuration you'll need.
You can also set up continuous integration on GitHub to check not only that
Sphinx can build the documentation in an isolated environment (i.e., with ``tox -e docs-test``)
but also that [ReadTheDocs can build it too](https://docs.readthedocs.io/en/stable/pull-requests.html).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/cinemol/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file
   configured to avoid the need for manual input at this step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion -- minor` after.

</details>