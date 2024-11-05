#!/usr/bin/env python3
import sys
from rdkit import Chem

def smiles2smart(smiles):
    m = Chem.MolFromSmiles(smiles)
    sma = Chem.MolToSmarts(m, isomericSmiles=True).replace('#6', 'C').replace('#8', 'O').replace('#15', 'P').replace('#7', 'N')
    print(sma)

if __name__ == "__main__":
    smiles2smart(sys.argv[1])