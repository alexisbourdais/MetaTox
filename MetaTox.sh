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
                 ===============
                ||             ||
                ||   MetaTox   ||
                ||             ||
                 ===============

 ==================================================
||                                                 ||
||   https://github.com/alexisbourdais/MetaTox/    ||
||                                                 ||
 ==================================================

####################
### Requirements ###
####################

- Singularity:
    sudo apt-get install -y singularity-container

- Conda :
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh

- Some package :
    sudo apt install -y zenity bc gawk dos2unix

- If error when creating metatrans conda environment, run before using MetaTox : conda config --set channel_priority flexible

#####################
###   Input file  ###
#####################

-> 1st column : Molecule ID/Name
-> 2nd column : SMILE
-> Separator  : comma ','

Example:

Nicotine,CN1CCC[C@H]1c2cccnc2


##################
###    Usage   ### 
##################

* With zenity : ./MetaTox.sh

* To skip zenity : ./MetaTox.sh --input input_file
    
    REQUIRED parameter

        -i|--input   

    OPTIONAL parameter

        -t|--type   Type of biotransformation to use with Biotransformer3:
                       [allHuman]  : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                        ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                        cyp450     : CYP450 metabolism prediction 
                        phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                        hgut       : Human gut microbial
                        superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                        envimicro  : Environmental microbial

        -n|--nstep  The number of steps for the prediction [default=2]

        -c|--cmode  CYP450 prediction Mode: 
                        1  = CypReact+BioTransformer rules
                        2  = CyProduct only
                       [3] = CypReact+BioTransformer rules+CyProducts
                    
        -1|--phase1 Number of reaction cycles Phase 1 [defaut=1]
        -2|--phase2 Number of reaction cycles Phase 2 [defaut=1]

"""
}

while [ $# -gt 0 ] ; do
    key="$1"
    case $key in
        -i|--input)
            input="$2"
            shift
            shift
            ;;
        -t|--type)
            type="$2"
            shift
            shift
            ;;
        -n|--nstep)
            nstep="$2"
            shift
            shift
            ;;
        -c|--cmode)
            cmode="$2"
            shift
            shift
            ;;
        -1|--phase1)
            phase1="$2"
            shift
            shift
            ;;
        -2|--phase2)
            phase2="$2"
            shift
            shift
            ;;
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

if conda info --envs | grep -q rdkit; then 
    echo "Conda environment "rdkit" already exists"
else 
    conda env create --name rdkit --file ${DirCondaEnv}rdkit_environment.yml
fi

if conda info --envs | grep -q metapredictor; then 
    echo "Conda environment "metapredictor" already exists"
else 
    conda env create --name metapredictor --file ${DirCondaEnv}metapred_environment.yml
fi

##############
### Zenity ###
##############

if [ -z $input ]; then

    zenity --info --text "
    In silico prediction by :
    - Bio-transformer3
    - Sygma
    - MetaTrans 
    - Meta-Predictor

    https://github.com/alexisbourdais/MetaTox/
    "

    input=$(zenity --file-selection --title="Select input file with each line = molecule,smile")
    option_meta=$(zenity --forms --title="MetaTrans options" --text="Directly Validate to apply default values" --add-entry="Min length (SMILE) [defaut=5]" --add-entry="Max length (SMILE) [defaut=120]" --add-entry="Top results [defaut=5 : top 10]" --separator=",")
        min=$(echo "$option_meta" | cut -d "," -f1)
        max=$(echo "$option_meta" | cut -d "," -f2)
        beam=$(echo "$option_meta" | cut -d "," -f3)

    option_sygma=$(zenity --forms --title="Sygma options" --text="Directly Validate to apply default values" --add-entry="Number of reaction cycles Phase 1 [defaut=1]" --add-entry="Number of reaction cycles Phase 2 [defaut=1]" --separator=",")
        phase1=$(echo "$option_meta" | cut -d "," -f1)
        phase2=$(echo "$option_meta" | cut -d "," -f2)

    type=$(zenity --list --title="Biotransformer 3 model" --text="Choose the type of biotransformation to use with Biotransformer3" --column="Type" --column="Description" \
    allHuman "Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step" ecbased "Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism" cyp450 "CYP450 metabolism prediction" phaseII "Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others" hgut "Human gut microbial" superbio "Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)" envimicro "Environmental microbial" )
    
    option_biotrans=$(zenity --forms --title="Biotransformer options" --text="Directly Validate to apply default values" --add-entry="The number of steps for the prediction [default=2]" --add-entry="CYP450 prediction Mode: 1=CypReact+BioTransformer rules; 2=CyProduct only; 3=CypReact+BioTransformer rules+CyProducts [Default=3]" --separator=",")
        nstep=$(echo "$option_biotrans" | cut -d "," -f1)
        cmode=$(echo "$option_biotrans" | cut -d "," -f2)
fi

###############
###  Input  ###
###############

file $input | grep CRLF && sudo dos2unix $input

declare -a tab_molecule
declare -a tab_smiles

while read a
do
    Molecule=$(echo $a | cut -d\, -f1)
    Smiles=$(echo $a | cut -d\, -f2)
    tab_smiles[${#tab_smiles[@]}]=${Smiles}
    tab_molecule[${#tab_molecule[@]}]=${Molecule}
done < $input

path_input="$(realpath ${input})"

##### Bio-Transformer option
###Default Mode
if [ -z $type ]; then
	type="allHuman"
fi
###Default nsteps
if [ -z $nstep ] || [[ ! $nstep = +([0-9]) ]]; then
	nstep="2"
fi
###Default cmode
if [ -z $cmode ] || [[ ! $cmode = +([0-3]) ]]; then
	cmode="3"
fi

#####Meta-Trans option
###Beam par défaut
if [ -z $BEAM ] || [[ ! $BEAM = +([0-9]) ]]; then
	BEAM=5
fi

###Min par défaut
if [ -z $MIN ] || [[ ! $MIN = +([0-9]) ]]; then
	MIN=5 
fi

###Max par défaut
if [ -z $MAX ] || [[ ! $MAX = +([0-9]) ]]; then
	MAX=120
fi

#####Sygma option
###phase1 par défaut
if [ -z $phase1 ] || [[ ! $phase1 = +([0-9]) ]]; then
	phase1=1
fi

###phase2 par défaut
if [ -z $phase2 ] || [[ ! $phase2 = +([0-9]) ]]; then
	phase2=1 
fi

###############
### Program ###
###############

for indice in ${!tab_molecule[@]}
do
    sygma_function --script ${DirScripts} --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput} -1 ${phase1} -2 ${phase2}
    biotransformer_function --script ${DirScripts} --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput} --biodir ${DirBiotrans} --type ${type} --nstep ${nstep} --cmode ${cmode}
done

metatrans_function --script ${DirScripts} --input "${path_input}" --outdir ${DirOutput} --metadir ${DirMetatrans} --min ${MIN} --max ${MAX} --beam ${BEAM}
metapredictor_function --script ${DirScripts} --input "${path_input}" --outdir ${DirOutput} --metadir ${DirMetapred}
    
    echo "
    Recording results in : ${DirOutput}

    Execution completed !
    "