#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def smiletoformula(smiles):
    mol = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(mol)
    print(formula)

if __name__ == "__main__":
    smiletoformula(sys.argv[1])