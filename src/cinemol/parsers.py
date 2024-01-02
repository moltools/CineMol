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
    
    atom_index = 0
    for line in lines[4:]:
                
        if line.strip() == "M  END":
            break

        parts = line.split()

        # Parse atom line.
        if len(parts) == 16:
            atom_index += 1
            x, y, z, atom_symbol = parts[0], parts[1], parts[2], parts[3]
            x, y, z = float(x), float(y), float(z)
            atoms.append(Atom(atom_index, atom_symbol, (x, y, z)))
        
        elif len(parts) == 7:
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