#!/usr/bin/env python3

import sys
from molmass import Formula

def mass_calcul(formule):
    Mass_hydro=1.007825
    Mass=Formula(formule).isotope.mass
    Mass_tot=Mass+Mass_hydro
    print(Mass_tot)

if __name__ == "__main__":
    mass_calcul(sys.argv[1])