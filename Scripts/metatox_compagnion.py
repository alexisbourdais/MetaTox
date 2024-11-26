#!/usr/bin/env python

###############
### Modules ###
###############
import argparse
import csv
import rdkit
from rdkit import Chem
import re
import os
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from molmass import Formula
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

#######################
### Argument parser ###
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-b", "--biotrans", required=True, help="biotransformer3 output")
ap.add_argument("-s", "--sygma", required=True, help="sygma output")
ap.add_argument("-mp", "--metapred", required=False, help="metapred output")
ap.add_argument("-mt", "--metatrans", required=False, help="metatrans output")
ap.add_argument("-o", "--output", required=True, help="output compilation")
ap.add_argument("-f", "--figure", required=True, help="List figure")
ap.add_argument("-d", "--dirfig", required=True, help="DirOutput figure")
args = vars(ap.parse_args())

#################
### Functions ###
#################
def smiles2smart(smiles):
    m = Chem.MolFromSmiles(smiles)
    sma = Chem.MolToSmarts(m, isomericSmiles=True).replace('#6', 'C').replace('#8', 'O').replace('#15', 'P').replace('#7', 'N').replace('#9', 'F')
    return(sma)

def smart2smile(smart):
    mol = Chem.rdmolfiles.MolFromSmarts(smart)
    smi = Chem.rdmolfiles.MolToSmiles(mol)
    return(smi)

def smiletoformula(smiles):
    mol = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(mol)
    return(formula)        

def mass_calcul(formule):
    Mass_hydro=1.007825
    Mass=Formula(formule).isotope.mass
    Mass_tot=Mass+Mass_hydro
    return(Mass_tot)

#The following functions are adapted from https://github.com/MunibaFaiza/cheminformatics/tree/main
def converter(file_name):
    sdf_file = Chem.SDMolSupplier(file_name)
    list_smiles=[]
    for mol in sdf_file:
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            list_smiles.append(smiles)
    return list_smiles

def smitostr(smile_txt, outDir):
    input=open(smile_txt, "r")
    for line in input:
        code_smile=line.split(",")
        smi=code_smile[1]
        molecule = Chem.MolFromSmiles(smi)
        fig = Draw.MolToMPL(molecule)
        plt.title(f"{code_smile[0]}")
        fig.savefig(f"{outDir}{code_smile[0]}.jpeg", bbox_inches='tight')

#################
### Variables ###
#################
smiles_list=[]
formulebrut_dic={}
mass_dic={}
biotrans_dic={}
biotrans_score_dic={}
biotrans_pathway_dic={}
biotrans_enzyme_dic={}
biotrans_system_dic={}
biotrans_precursor_smiles_dic={}
biotrans_precursor_formule_dic={}
sygma_dic={}
sygma_pathway_dic={}
sygma_score_dic={}
metatrans_dic={}
metapred_dic={}

#######################
### BioTransformer3 ###
#######################
biotrans = csv.reader(open(args['biotrans'], "r"), delimiter=';')

for line in biotrans:
    if "SMILES" in line:
        continue
    else:
        smiles = line[2]
        smart=smiles2smart(smiles)
        new_smiles=smart2smile(smart)
        score = line[7]
        pathway = line[13]
        enzyme = line[15]
        system = line[16]

        try:
            formulebrute=smiletoformula(new_smiles)
            mass=mass_calcul(formulebrute)
        except:
            formulebrute="NA"
            mass="NA"

        smiles_precursor = line[18]
        smart_precursor=smiles2smart(smiles_precursor)
        new_smiles_precursor=smart2smile(smart_precursor)

        try:
            formulebrute_precursor=smiletoformula(new_smiles_precursor)
            #mass_precursor=mass_calcul(formulebrute)

        except:
            formulebrute_precursor="NA"
            #mass_precursor="NA"

        smiles_list.append(new_smiles)
        formulebrut_dic.setdefault(new_smiles, formulebrute)
        mass_dic.setdefault(new_smiles, mass)
        biotrans_dic.setdefault(new_smiles, "+")
        biotrans_score_dic.setdefault(new_smiles, score)
        biotrans_pathway_dic.setdefault(new_smiles, pathway)
        biotrans_enzyme_dic.setdefault(new_smiles, enzyme)
        biotrans_system_dic.setdefault(new_smiles, system)
        biotrans_precursor_smiles_dic.setdefault(new_smiles, new_smiles_precursor)
        biotrans_precursor_formule_dic.setdefault(new_smiles, formulebrute_precursor)

#############
### SygMa ###
#############
sdf_smiles=args['sygma']

list_smiles_sygma=converter(sdf_smiles)

list_pathway=[]
list_score=[]
line_pre=""

#Get score and pathway from sdf file
for line in open(sdf_smiles, 'r'):
    if line_pre == "pathway":
        pathway=line.replace(",", "~").strip()
        line_pre="pathway_debut"
        continue

    if line_pre == "pathway_debut":
        pathway_suite=line.replace(",", "~").replace(";", "").strip()
        pathway_tot="".join([pathway, pathway_suite])
        list_pathway.append(pathway_tot)
        line_pre=""

    if line_pre == "score":
        score=line.strip()
        list_score.append(score)
        line_pre=""

    if re.match("^>  <Pathway>", line):
        line_pre="pathway"

    if re.match("^>  <Score>", line):
        line_pre="score"

for i in range(len(list_smiles_sygma)):
    smiles=list_smiles_sygma[i]
    score=list_score[i]
    pathway=list_pathway[i]

    if pathway == "parent":
        pass
    else:
        smart=smiles2smart(smiles)
        new_smiles=smart2smile(smart)
        sygma_dic.setdefault(new_smiles, "+")
        sygma_pathway_dic.setdefault(new_smiles, pathway)
        sygma_score_dic.setdefault(new_smiles, score)

        if new_smiles in smiles_list:
            pass
        else:
            try:
                formulebrute=smiletoformula(new_smiles)
                mass=mass_calcul(formulebrute)
            except:
                formulebrute="NA"
                mass="NA"

            smiles_list.append(new_smiles)
            formulebrut_dic.setdefault(new_smiles, formulebrute)
            mass_dic.setdefault(new_smiles, mass)

##################
### Meta-Trans ###
##################
if os.path.isfile(args['metatrans']):
    metatrans_file=open(args['metatrans'], "r")

    for line in metatrans_file:
        smiles = line
        smart=smiles2smart(smiles)
        new_smiles=smart2smile(smart)
        metatrans_dic.setdefault(new_smiles, "+")

        if new_smiles in smiles_list:
            pass
        else:
            try:
                formulebrute=smiletoformula(new_smiles)
                mass=mass_calcul(formulebrute)
            except:
                formulebrute="NA"
                mass="NA"

            smiles_list.append(new_smiles)
            formulebrut_dic.setdefault(new_smiles, formulebrute)
            mass_dic.setdefault(new_smiles, mass)

#################
### Meta-Pred ###
#################
if os.path.isfile(args['metapred']):
    metapred_file=open(args['metapred'], "r")

    for line in metapred_file:
        smiles = line
        smart=smiles2smart(smiles)
        new_smiles=smart2smile(smart)
        metapred_dic.setdefault(new_smiles, "+")

        if new_smiles in smiles_list:
            pass
        else:
            try:
                formulebrute=smiletoformula(new_smiles)
                mass=mass_calcul(formulebrute)
            except:
                formulebrute="NA"
                mass="NA"
        
            smiles_list.append(new_smiles)
            formulebrut_dic.setdefault(new_smiles, formulebrute)
            mass_dic.setdefault(new_smiles, mass)

###################
### Compilation ###
###################
results_file=open(args['output'], "a")
figures_file=args['figure']
dir_figures=args['dirfig']

#Entete
print("FormuleBrute\tMasse(+H)\tSmiles\tSygma\tBioTransformer3\tMetaTrans\tMetaPredictor\tSygma_pathway\tBioTrans_pathway\tSygma_score\tBioTrans_score\tBioTrans_precursor\tBioTrans_precursor\tBioTrans_enzyme\tBioTrans_system\tFigure", file=results_file)

#Count metabolites
nbmolecule=0

for smiles in smiles_list:
    formulebrute=formulebrut_dic.get(smiles)
    mass=mass_dic.get(smiles)
    biotrans=biotrans_dic.get(smiles)
    biotrans_score=biotrans_score_dic.get(smiles)
    biotrans_pathway=biotrans_pathway_dic.get(smiles)
    biotrans_enzyme=biotrans_enzyme_dic.get(smiles)
    biotrans_system=biotrans_system_dic.get(smiles)
    biotrans_prec_smiles=biotrans_precursor_smiles_dic.get(smiles)
    biotrans_prec_for=biotrans_precursor_formule_dic.get(smiles)
    sygma=sygma_dic.get(smiles)
    sygma_pathway=sygma_pathway_dic.get(smiles)
    sygma_score=sygma_score_dic.get(smiles)
    metatrans=metatrans_dic.get(smiles)
    metapred=metapred_dic.get(smiles)

    #Files to create figures
    if formulebrute == "NA":
        figure="NA"
    else:
        nbmolecule+=1
        figure=f"Figure_{nbmolecule}"
        print(f"Molecule_{nbmolecule},{smiles}", file=open(figures_file, 'a'))
    
    print(f"{formulebrute}\t{mass}\t{smiles}\t{sygma}\t{biotrans}\t{metatrans}\t{metapred}\t{sygma_pathway}\t{biotrans_pathway}\t{sygma_score}\t{biotrans_score}\t{biotrans_prec_for}\t{biotrans_prec_smiles}\t{biotrans_enzyme}\t{biotrans_system}\t{figure}".replace("None", ""), file=results_file)

#Figures creation
smitostr(figures_file, dir_figures)