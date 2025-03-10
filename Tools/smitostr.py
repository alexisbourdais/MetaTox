#!/usr/bin/env python3
import argparse
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--infile", required=True, help="input file")
args = vars(ap.parse_args())

with open(args['infile'], "r") as smile_txt:
        for line in smile_txt:
            code_smile=line.split(",")
            smi=code_smile[1]
            molecule = Chem.MolFromSmiles(smi)
            fig = Draw.MolToMPL(molecule)
            plt.title(f"{code_smile[0]}")
            fig.savefig(f"{code_smile[0]}.jpeg", bbox_inches='tight')