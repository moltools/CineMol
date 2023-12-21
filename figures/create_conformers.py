#!/usr/bin/env python3
"""
https://www.npatlas.org/download
"""
import typing as ty
import argparse 
import errno 
import functools
import os 
import signal
from rdkit import Chem 
from rdkit.Chem import AllChem

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-npatlas_tsv", type=str, required=True)
    parser.add_argument("-out", type=str, required=True)
    return parser.parse_args()

class TimeoutError(Exception):
    pass 

def timeout(seconds: int = 5, error_message: str = os.strerror(errno.ETIME)) -> ty.Callable:
    """
    Decorator to raise a TimeoutError when runtime exceeds the specified time.

    Parameters
    ----------
    seconds : int, optional
        Timeout in seconds, by default 5.
    error_message : str, optional
        Error message to be raised, by default os.strerror(errno.ETIME).
    
    Note: Timer only works on UNIX systems; also not thread-safe.
    """
    def decorator(func: ty.Callable) -> ty.Callable:

        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try: 
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result
        
        return wrapper
    
    return decorator

@timeout(seconds=3)
def embed_conformer(mol: Chem.Mol) -> ty.Optional[Chem.Mol]:
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xf00d, maxAttempts=1000)
    AllChem.MMFFOptimizeMolecule(mol) # MMFF94
    try: _ = mol.GetConformer()
    except: return None
    return mol

def main() -> None:
    args = cli()
    parsed_mols = 0
    writer = Chem.SDWriter(args.out)
    with open(args.npatlas_tsv, "r") as f:
        f.readline() # skip header
        for i, line in enumerate(f):
            smiles = line.strip().split("\t")[10]
            try:
                mol = Chem.MolFromSmiles(smiles)
                heavy_atom_count = mol.GetNumHeavyAtoms()
                if heavy_atom_count > 50: continue
                mol = embed_conformer(mol)
                if mol is not None:
                    writer.write(mol)
                    parsed_mols += 1
            except:
                continue
            print(f"{i}".zfill(10) + " " + f"{parsed_mols}".zfill(10), end="\r")
    writer.close()
    exit(0)

if __name__ == "__main__":
    main()