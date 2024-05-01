"""
Description:    Measure performance of CineMol's algorithm.
Usage:          python3 measure_performance.py -i path/to/sdf/file.sdf -o path/to/out/file.tsv
"""

import argparse
import time
import typing as ty

from rdkit import Chem

from cinemol.api import Atom, Bond, Look, Style, draw_molecule


def cli() -> argparse.Namespace:
    """
    Command line interface for this script.

    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to input SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Path to output TSV file.")
    return parser.parse_args()


def parse_mol(mol: Chem.Mol) -> ty.Tuple[ty.List[Atom], ty.List[Bond]]:
    """
    Parse a molecule into a list of atoms and bonds. Exclude hydrogens.

    :param Chem.Mol mol: Molecule to parse.
    :return: Tuple of lists of atoms and bonds.
    :rtype: ty.Tuple[ty.List[Atom], ty.List[Bond]]
    """
    pos = mol.GetConformer().GetPositions()

    atoms, bonds = [], []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        atoms.append(Atom(atom.GetIdx(), atom.GetSymbol(), pos[atom.GetIdx()]))

    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol() == "H" or bond.GetEndAtom().GetSymbol() == "H":
            continue
        start_index, end_index = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds.append(Bond(start_index, end_index, int(bond.GetBondTypeAsDouble())))

    return atoms, bonds


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    out_file = open(args.o, "w")
    out_file.write("index\tnum_heavy_atoms\tnum_bonds\tdepiction\tlook\ttime\tsize\n")

    # Parse mols from SDF file with RDKit.
    mols = [mol for mol in Chem.SDMolSupplier(args.i)]

    # Parse mols, draw with different depictions and styles, and report runtime and file size.
    for i, mol in enumerate(mols):
        num_heavy_atoms = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() != "H"])
        num_bonds = len(
            [
                bond
                for bond in mol.GetBonds()
                if bond.GetBeginAtom().GetSymbol() != "H" and bond.GetEndAtom().GetSymbol() != "H"
            ]
        )
        atoms, bonds = parse_mol(mol)

        styles = [Style.SpaceFilling, Style.BallAndStick, Style.Tube, Style.Wireframe]
        looks = [Look.Cartoon, Look.Glossy]

        for look in looks:
            for style in styles:
                t0 = time.time()
                svg = draw_molecule(atoms=atoms, bonds=bonds, style=style, look=look, resolution=50)
                svg_str = svg.to_svg()

                runtime_ms = 1000 * (time.time() - t0)
                file_size_kb = len(svg_str) / 1000

                out_file.write(
                    f"{i}\t{num_heavy_atoms}\t{num_bonds}\t{style.name}\t{look.name}\t{runtime_ms}\t{file_size_kb}\n"
                )
                out_file.flush()

        padding = len(str(len(mols)))
        print(f"{i}".zfill(padding), end="\r")

    out_file.close()

    exit(0)


if __name__ == "__main__":
    main()
