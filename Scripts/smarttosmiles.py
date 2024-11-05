#!/usr/bin/env python3
import sys
from rdkit import Chem

def smart2smile(smart):

    mol = Chem.rdmolfiles.MolFromSmarts(smart)
    smi = Chem.rdmolfiles.MolToSmiles(mol)
    print(smi)

if __name__ == "__main__":
    smart2smile(sys.argv[1])