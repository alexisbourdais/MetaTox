#!  /usr/bin/bash

#####################
### Help Function ###
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
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh

- Some packages :
    sudo apt install -y zenity gawk dos2unix csvkit

- If error when creating metatrans conda environment, run before using MetaTox : conda config --set channel_priority flexible

#####################
###   Input file  ###
#####################

-> 1st column : Molecule ID/Name
-> 2nd column : SMILES
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

        -m|--meta      To activate metaTrans and meta-Predictor [No]

        -t|--type       Type of biotransformation to use with BioTransformer3:
                            [allHuman] : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                            ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                            cyp450     : CYP450 metabolism prediction 
                            phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                            hgut       : Human gut microbial
                            superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                            envimicro  : Environmental microbial

        -n|--nstep      The number of steps for the prediction by BioTransformer [default=1]

        -c|--cmode      CYP450 prediction Mode uses by BioTransformer: 
                            1  = CypReact+BioTransformer rules
                            2  = CyProduct only
                           [3] = CypReact+BioTransformer rules+CyProducts
                    
        -1|--phase1     Number of reaction cycles Phase 1 by SygMa [defaut=1]
        -2|--phase2     Number of reaction cycles Phase 2 by SygMa [defaut=1]

        -p|--tmp        To keep intermediate files [No]

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
        -m|--meta)
            meta_activate=true
            shift
            ;;
        -p|--tmp)
            keep_tmp=true
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

tmp="${work_dir}/tmp/"
if test -d $tmp; then
  rm -r $tmp
fi
mkdir $tmp

#For download version of biotransformer3
#DirBiotrans="${work_dir}/biotransformer3.0jar/"
DirMetatrans="${work_dir}/MetaTrans/"
DirMetapred="${work_dir}/Meta-Predictor/"
DirCondaEnv="${work_dir}/CondaEnv/"
DirScripts="${work_dir}/Scripts/"

Script_Metatox_companion="${DirScripts}metatox_compagnion.py"
Script_SmitoStr="${DirScripts}smitostr.py"

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

    input=$(zenity --file-selection --title="Select input file with each line = molecule,smiles")
    option_meta=$(zenity --forms --title="MetaTrans options" --text="Directly Validate to apply default values" --add-entry="Min length (SMILES) [defaut=5]" --add-entry="Max length (SMILES) [defaut=120]" --add-entry="Top results [defaut=5 : top 10]" --separator=",")
        min=$(echo "$option_meta" | cut -d "," -f1)
        max=$(echo "$option_meta" | cut -d "," -f2)
        beam=$(echo "$option_meta" | cut -d "," -f3)

    option_sygma=$(zenity --forms --title="Sygma options" --text="Directly Validate to apply default values" --add-entry="Number of reaction cycles Phase 1 [defaut=1]" --add-entry="Number of reaction cycles Phase 2 [defaut=1]" --separator=",")
        phase1=$(echo "$option_sygma" | cut -d "," -f1)
        phase2=$(echo "$option_sygma" | cut -d "," -f2)

    type=$(zenity --list --title="Biotransformer 3 model" --text="Choose the type of biotransformation to use with Biotransformer3" --column="Type" --column="Description" \
    allHuman "Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step" ecbased "Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism" cyp450 "CYP450 metabolism prediction" phaseII "Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others" hgut "Human gut microbial" superbio "Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)" envimicro "Environmental microbial" )
    
    option_biotrans=$(zenity --forms --title="Biotransformer options" --text="Directly Validate to apply default values" --add-entry="The number of steps for the prediction [default=1]" --add-entry="CYP450 prediction Mode: 1=CypReact+BioTransformer rules; 2=CyProduct only; 3=CypReact+BioTransformer rules+CyProducts [Default=3]" --separator=",")
        nstep=$(echo "$option_biotrans" | cut -d "," -f1)
        cmode=$(echo "$option_biotrans" | cut -d "," -f2)

    zenity --question --title="Meta-Predictor Activation" --text="Do you want to activate MetaTrans and Meta-Predictor ?"
    metapred_activate_answer=$?
    if [ $metapred_activate_answer -eq 0 ]; then
	    meta_activate=true
    else [ $metapred_activate_answer -eq 1 ]
    	meta_activate=false
    fi
fi

########################
### Defaults Options ###
########################

##### Bio-Transformer options
###Default Mode
if [ -z $type ]; then
	type="allHuman"
fi

###Default nsteps
if [ -z $nstep ] || [[ ! $nstep = +([0-9]) ]]; then
	nstep="1"
fi

###Default cmode
if [ -z $cmode ] || [[ ! $cmode = +([0-3]) ]]; then
	cmode="3"
fi

##### Meta-Trans options
###Default Beam
if [ -z $BEAM ] || [[ ! $BEAM = +([0-9]) ]]; then
	BEAM=5
fi

###Default min
if [ -z $MIN ] || [[ ! $MIN = +([0-9]) ]]; then
	MIN=5 
fi

###Default Max
if [ -z $MAX ] || [[ ! $MAX = +([0-9]) ]]; then
	MAX=120
fi

##### Sygma options
###Default phase1
if [ -z $phase1 ] || [[ ! $phase1 = +([0-9]) ]]; then
	phase1=1
fi

###Default phase2
if [ -z $phase2 ] || [[ ! $phase2 = +([0-9]) ]]; then
	phase2=1 
fi

##### Meta-Predictor option
###Activation
if [ -z $meta_activate ]; then
	meta_activate=false
fi

###Keep intermediate files
if [ -z $keep_tmp ]; then
	keep_tmp=false
fi

echo "
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

    - BioTransformer bio-reaction : $type
    - BioTransformer cycle number : $nstep
    - Biotransformer CYP450 prediction Mode : $cmode

    - SygMa Phase 1 cycle number : $phase1
    - SygMa Phase 2 cycle number : $phase2

    - MetaTrans & Meta-Predictor activation : $meta_activate
"

#########################
### Environment conda ###
#########################

eval "$(conda shell.bash hook)"

if conda info --envs | grep -q rdkit; then 
    :
else
    echo "Installation of "rdkit" environment :"
    conda env create --name rdkit --file ${DirCondaEnv}rdkit_environment.yml
fi

if $meta_activate; then
    if conda info --envs | grep -q metapredictor; then 
        :
    else 
        echo "Installation of "metapredictor" environment :"
        conda env create --name metapredictor --file ${DirCondaEnv}metapred_environment.yml
    fi

    if conda info --envs | grep -q metatrans; then 
        :
    else
        echo "Installation of "metatrans" environment :"
        conda env create --name metatrans --file ${DirCondaEnv}metatrans_environment.yml
    fi
fi

###############
###  Input  ###
###############

file $input | grep CRLF && dos2unix $input

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

if ${meta_activate}; then

    ##################
    ### META-TRANS ###
    ##################

    conda activate metatrans

    cd $DirMetatrans

    #Variables from MetaTrans
    outfile=processed_data.txt
    results=Prediction_MetaTrans.csv
    images=./Figures/*
    STORE=predictions/  #directory for output process
    mkdir $STORE

    echo "
     ==================================================
    ||                                                ||
    ||   Process of MetaTrans step 1 : Prepara data   ||
    ||                                                ||
     ==================================================
    "

    python prepare_input_file.py -input_file "${path_input}" -output_file ${outfile}

    echo "
     ==================================================
    ||                                                ||
    ||   Process of MetaTrans step 2 : Translate      ||
    ||                                                ||
     ==================================================
    "

    for model_id in {1,2,3,4,5,6}
    do
        MODEL_FILE='models/model_'$model_id'.pt'
        OUT_NAME='model'$model_id'_beam'$BEAM'.txt'
        OUT_FILE=$STORE$OUT_NAME
        python translate.py -model $MODEL_FILE -src $outfile -output $OUT_FILE -n_best $BEAM -beam_size $BEAM -verbose -min_length $MIN -max_length $MAX > "${tmp}MetaTrans_log.txt" 2>&1
    done

    echo "
     ========================================================
    ||                                                      ||
    ||   Process of MetaTrans step 3 : Get prediction       ||
    ||                                                      ||
     ========================================================
    "

    python process_predictions.py \
    -input_file "${path_input}" \
    -output_file ${results} \
    -beam_size ${BEAM} \
    -visualise_molecules True \
    >> "${tmp}MetaTrans_log.txt" 2>&1

    conda deactivate

    #Moving results files and deleting intermediate files
    rm -r ${outfile} ${STORE} ${images}
    mv ${results} ${tmp}
    cd ${tmp}

    #Separation of results into separate files for each molecule
    pat='^Molecule'

    while read a
    do
        if [[ $a =~ $pat ]]; then
            :

        else
            MoleculeID=$(echo $a | cut -d\, -f1)
            #SmileParent=$(echo $a | cut -d\, -f2)
            SmileMetabo=$(echo $a | cut -d\, -f3)
            read -ra tab_metabo <<< "$SmileMetabo"
            for i in ${tab_metabo[@]}; do
                echo "${i}" >> "${MoleculeID}_Metatrans.csv"
            done
        fi
    done < "Prediction_MetaTrans.csv"

    ######################
    ### META-PREDICTOR ###
    ######################

    echo "
     =================================
    ||                               ||
    ||   Process of MetaPredictor    ||
    ||                               ||
     =================================
    "

    conda activate metapredictor

    cd ${DirMetapred}

    python prepare_input_file.py -input_file ${path_input} -output_file processed_data.txt
    bash predict-top15.sh processed_data.txt ./prediction ${path_input} > "${tmp}Metapredictor_log" 2>&1

    conda deactivate

    ### Mv results and remove temp files
    mv prediction/predict.csv ${tmp}Prediction_metapred.csv
    rm -r prediction/* Figures/* processed_data.txt

    #Separation of results into separate files for each molecule
    cd ${tmp}
    pat='^Name'

    while read a
    do
        if [[ $a =~ $pat ]]; then
            :

        else
            MoleculeID=$(echo $a | cut -d\, -f1)
            SmileMetabo=$(echo $a | cut -d\, -f3)
            read -ra tab_metabo <<< "$SmileMetabo"
            for i in ${tab_metabo[@]}; do
                echo "${i}" >> "${MoleculeID}_metapred.csv"
            done
        fi
    done < "Prediction_metapred.csv"
fi

for indice in ${!tab_molecule[@]}
do
    results_file="${DirOutput}${tab_molecule[${indice}]}_CompileResults.tsv"
    results_figure="${DirOutput}${tab_molecule[${indice}]}_figures/"
    mkdir ${results_figure}

    #########################
    ### BIOTRANSFORMERS 3 ###
    #########################

    echo "
     ===============================================================
    ||                                                                                           
    ||   Process of ${tab_molecule[${indice}]} by Biotransformer3
    ||                                                                                            
     ===============================================================
    "
    #Download version
    #java -jar BioTransformer3.0_20230525.jar -b ${type} -k pred -cm 3 -ismi "${smile}" -ocsv "${DirOutput}${molecule}_Biotransformer3.csv"

    #Singularity version
    singularity exec https://depot.galaxyproject.org/singularity/biotransformer%3A3.0.20230403--hdfd78af_0 biotransformer \
    -b "${type}" \
    -k "pred" \
    -cm 3 \
    -s "${nstep}" \
    -ismi "${tab_smiles[${indice}]}" \
    -ocsv "${tmp}${tab_molecule[${indice}]}_Biotransformer3.csv" 2>&1 | tee -a "${tmp}${tab_molecule[${indice}]}_Biotransformer3_log.txt"

    #Changement csv format
    csvformat -D ";" "${tmp}${tab_molecule[${indice}]}_Biotransformer3.csv" | gawk -v RS='"' 'NR % 2 == 0 { gsub(/\n/, "") } { printf("%s%s", $0, RT) }' > "${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut.csv"

    ##################
    ###    SygMa   ###
    ##################

    echo "
     ====================================================
    ||                                         
    ||   Process of ${tab_molecule[${indice}]} by Sygma
    ||                                             
     ====================================================
    "

    singularity run docker://3dechem/sygma ${tab_smiles[${indice}]} -1 $phase1 -2 $phase2 >> "${tmp}${tab_molecule[${indice}]}_Sygma.sdf"

    ###################
    ### Compilation ###
    ###################

    echo "     
     =================================
    ||                               ||
    ||      Results Compilation      ||
    ||                               ||
     =================================
    "

    conda activate rdkit

    python3 $Script_Metatox_companion \
    --biotrans "${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut.csv" \
    --sygma "${tmp}${tab_molecule[${indice}]}_Sygma.sdf" \
    --metapred "${tmp}${tab_molecule[${indice}]}_metapred.csv" \
    --metatrans "${tmp}${tab_molecule[${indice}]}_Metatrans.csv" \
    --output "${results_file}" \
    --figure "${tmp}${tab_molecule[${indice}]}_ListeSmile.txt" \
    --dirfig "${results_figure}" \
    > "${tmp}${tab_molecule[${indice}]}_Compagnion_log.txt" 2>&1

    conda deactivate
done

if ${keep_tmp}; then
    :
else
    rm -r $tmp
fi

echo "
Recording results in : ${DirOutput}

Execution completed !
"