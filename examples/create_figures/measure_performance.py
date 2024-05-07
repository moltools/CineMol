"""
Description:    Measure performance of CineMol's algorithm.
Dependencies:   rdkit==2023.9.0; tqdm==4.66.4
Usage:          python3 measure_performance.py -i path/to/sdf/file.sdf -o path/to/out/file.tsv
"""

import argparse
import multiprocessing as mp
import time
import typing as ty

from rdkit import Chem
from tqdm import tqdm

from cinemol.api import Atom, Bond, Look, Style, draw_molecule


class Job:
    """Record of a job to draw a molecule with different depictions and styles."""

    def __init__(self, i: int, mol: Chem.Mol, style: Style, look: Look):
        """
        Initialize a job.

        :param int i: Index of the molecule.
        :param Chem.Mol mol: Molecule to draw.
        :param Style style: Style of the depiction.
        :param Look look: Look of the depiction.
        """
        self.i = i
        self.mol = mol
        self.style = style
        self.look = look


def cli() -> argparse.Namespace:
    """
    Command line interface for this script.

    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to input SDF file.")
    parser.add_argument("-o", type=str, required=True, help="Path to output TSV file.")
    parser.add_argument("-n", type=int, default=1, help="Number of threads to use.")
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


def run_job(job: Job) -> ty.Tuple[int, int, int, float, float]:
    """
    Run a job to draw a molecule with different depictions and styles.

    :param job: Job to run.
    :type job: Job
    :return: Tuple of index, number of heavy atoms, number of bonds, runtime in milliseconds, and file size in kilobytes.
    :rtype: ty.Tuple[int, int, int, float, float]
    """
    i, mol, style, look = job.i, job.mol, job.style, job.look

    num_heavy_atoms = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() != "H"])
    num_bonds = len(
        [
            bond
            for bond in mol.GetBonds()
            if bond.GetBeginAtom().GetSymbol() != "H" and bond.GetEndAtom().GetSymbol() != "H"
        ]
    )
    atoms, bonds = parse_mol(mol)

    t0 = time.time()
    svg = draw_molecule(
        atoms=atoms, bonds=bonds, style=style, look=look, resolution=50, scale=1000.0
    )
    svg_str = svg.to_svg()

    runtime_ms = 1000 * (time.time() - t0)
    file_size_kb = len(svg_str) / 1000

    return i, num_heavy_atoms, num_bonds, runtime_ms, file_size_kb, style, look


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    out_file = open(args.o, "w")
    out_file.write("index\tnum_heavy_atoms\tnum_bonds\tdepiction\tlook\ttime\tsize\n")

    # Parse mols from SDF file with RDKit.
    mols = [mol for mol in Chem.SDMolSupplier(args.i)]

    # Jobs.
    jobs = []

    # Determine number of threads to use.
    num_threads = min(args.n, mp.cpu_count())

    # Create jobs.
    for i, mol in enumerate(mols):
        styles = [Style.SPACEFILLING, Style.BALL_AND_STICK, Style.TUBE, Style.WIREFRAME]
        looks = [Look.CARTOON, Look.GLOSSY]

        for look in looks:
            for style in styles:
                jobs.append(Job(i, mol, style, look))

    # Run jobs multi-threaded.
    with mp.Pool(processes=num_threads) as pool:
        for i, result in tqdm(enumerate(pool.imap_unordered(run_job, jobs))):

            i, num_heavy_atoms, num_bonds, runtime_ms, file_size_kb, style, look = result

            out_file.write(
                f"{i}\t{num_heavy_atoms}\t{num_bonds}\t{style.name}\t{look.name}\t{runtime_ms}\t{file_size_kb}\n"
            )
            out_file.flush()

    out_file.close()

    exit(0)


if __name__ == "__main__":
    main()
