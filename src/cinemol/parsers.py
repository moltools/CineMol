import typing as ty 

from cinemol.chemistry import Atom, Bond 

def parse_sdf(src: str, include_hs: bool = True) -> ty.Tuple[ty.List[Atom], ty.List[Bond]]:
    """
    Parse first molecule from SDF file format.

    :param string src: SDF file content.
    :param bool include_hs: Include hydrogens.
    :return: Atoms and bonds.
    :rtype: ty.Tuple[ty.List[Atom], ty.List[Bond]]
    """
    atoms, bonds = [], []

    lines = src.split("$$$$\n")[0].split("\n")

    counts_line = lines[3]
    atom_count_str, bond_count_str, *_ = counts_line.strip().split()
    atom_count, bond_count = int(atom_count_str), int(bond_count_str)

    atom_lines = lines[4:4 + atom_count]
    bond_lines = lines[4 + atom_count:4 + atom_count + bond_count]

    atom_index = 0

    # Parse atom line.
    for atom_line in atom_lines:
        atom_index += 1
        parts = atom_line.strip().split()
        x, y, z, atom_symbol = parts[0], parts[1], parts[2], parts[3]
        x, y, z = float(x), float(y), float(z)
        atoms.append(Atom(atom_index, atom_symbol, (x, y, z)))
        
    # Parse bond line.
    for bond_line in bond_lines:
        parts = bond_line.strip().split()
        start_index, stop_index, bond_order = parts[0], parts[1], parts[2]
        bonds.append(Bond(int(start_index), int(stop_index), int(bond_order)))

    atom_map = {atom.index: atom for atom in atoms}
    if not include_hs:
        atoms = [atom for atom in atoms if atom.symbol != "H"]
        bonds = [
            bond for bond in bonds 
            if (
                atom_map[bond.start_index].symbol != "H" and 
                atom_map[bond.end_index].symbol != "H"
            )
        ]

    return atoms, bonds