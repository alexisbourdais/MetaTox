#!  /usr/bin/bash

source ${PWD}/Scripts/biotransformer3.sh
source ${PWD}/Scripts/metapredictor.sh
source ${PWD}/Scripts/metatrans.sh
source ${PWD}/Scripts/sygma.sh

#####################
### Fonction HELP ###
#####################

help_msg() { 
    printf """

Usage: ./$0

#####################
###   Input file  ###
#####################

-> 1ère colonne : ID/Nom molecule
-> 2eme colonne : code SMILE
-> Séparateur   : virgule ','

Exemple:

Nicotine,CN1CCC[C@H]1c2cccnc2

####################
### Requirements ###
####################

- Singularity pour l'outil Sygma (https://docs.docker.com/)
- Conda pour l'outil MetaTrans (https://docs.conda.io/projects/conda/en/stable/user-guide/index.html)

"""
}

while [ $# -gt 0 ] ; do
    key="$1"
    case $key in
        *)
            help_msg
            exit 0                         
    esac
done

######################
### Work Directory ###
######################

work_dir="${PWD}"

DirOutput="${work_dir}/Results_Prediction/"
mkdir $DirOutput

DirScripts="${work_dir}/Scripts/"
DirBiotrans="${work_dir}/biotransformer3.0jar/"
DirMetatrans="${work_dir}/MetaTrans/"
DirMetapred="${work_dir}/Meta-Predictor/"
DirCondaEnv="${work_dir}/CondaEnv/"

#########################
### Environment conda ###
#########################

eval "$(conda shell.bash hook)"
#conda init --all

if conda info --envs | grep -q metatrans; then 
    echo "Conda environment "metatrans" already exists"
else 
    conda env create --name metatrans --file ${DirCondaEnv}metatrans_environment.yml
fi

if conda info --envs | grep -q my-rdkit-env; then 
    echo "Conda environment "my-rdkit-env" already exists"
else 
    conda env create --name my-rdkit-env --file ${DirCondaEnv}rdkit_environment.yml
fi

### Add auto installation if possible
if conda info --envs | grep -q metapredictor; then 
    echo "Conda environment "metapredictor" already exists"
else 
    echo "Conda environment "metapredictor" does not exists !
        You need to install it manually (https://github.com/zhukeyun/Meta-Predictor)
        "
    exit 0
fi

##############
### Zenity ###
##############

zenity --info --text "
In silico prediction by :
- Bio-transformer3
- Sygma
- MetaTrans 
- Meta-Predictor
"

input_select=$(zenity --file-selection --title="Select input file with each line = molecule,smile")
option_meta=$(zenity --forms --title="MetaTrans options" --text="Directly Validate to apply default values" --add-entry="Min length (SMILE) [defaut=5]" --add-entry="Max length (SMILE) [defaut=120]" --add-entry="Top results [defaut=5 : top 10]" --separator=",")
min=$(echo "$option_meta" | cut -d "," -f1)
max=$(echo "$option_meta" | cut -d "," -f2)
beam=$(echo "$option_meta" | cut -d "," -f3)
type=$(zenity --list --title="Biotransformer 3 model" --text="Choose the model to use with Biotransformer 3" --column="Model" --column="Description" \
allHuman "Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step" ecbased "Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism" cyp450 "CYP450 metabolism prediction" phaseII "Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others" hgut "Human gut microbial" superbio "Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)" envimicro "Environmental microbial" )

###############
###  Input  ###
###############

file $input_select | grep CRLF && sudo dos2unix $input_select

declare -a tab_molecule
declare -a tab_smiles

while read a
do
    Molecule=$(echo $a | cut -d\, -f1)
    Smiles=$(echo $a | cut -d\, -f2)
    tab_smiles[${#tab_smiles[@]}]=${Smiles}
    tab_molecule[${#tab_molecule[@]}]=${Molecule}
done < $input_select

###############
### Program ###
###############

for indice in ${!tab_molecule[@]}
do
    sygma_function --script $DirScripts --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput}
    biotransformer_function --script $DirScripts --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput} --biodir $DirBiotrans --type $type
done

metatrans_function --script $DirScripts --input $input_select --outdir ${DirOutput} --metadir ${DirMetatrans} --min $min --max $max --beam $beam

if conda info --envs | grep -q metapredictor; then 
    metapredictor_function --script $DirScripts --input $input_select --outdir ${DirOutput} --metadir ${DirMetapred}
    echo "
    Enregistrement des résultats dans le dossier ${DirOutput}

    Fin de l'excecution du programme !
    "
else 
    echo "Conda environment "metapredictor" does not exists !
    You need to install it manually (https://github.com/zhukeyun/Meta-Predictor)

    Enregistrement des résultats dans le dossier ${DirOutput}

    Fin de l'excecution du programme !
    "

    exit 0
fi