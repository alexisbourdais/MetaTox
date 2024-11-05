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
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh

- Some packages :
    sudo apt install -y zenity bc gawk dos2unix

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

        -m|--metapred   To activate metapredictor [No]

        -t|--type       Type of biotransformation to use with BioTransformer3:
                            [allHuman]  : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                            ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                            cyp450     : CYP450 metabolism prediction 
                            phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                            hgut       : Human gut microbial
                            superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                            envimicro  : Environmental microbial

        -n|--nstep      The number of steps for the prediction by BioTransformers [default=2]

        -c|--cmode      CYP450 prediction Mode uses by BioTransformers: 
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
        -m|--metapred)
            metapred_activate=true
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

#####################
### Progress Tool ###
#####################

bar_size=40
bar_char_done="#"
bar_char_todo="-"
bar_percentage_scale=2

function show_progress {
    current="$1"
    total="$2"

    # calculate the progress in percentage 
    percent=$(bc <<< "scale=$bar_percentage_scale; 100 * $current / $total" )
    # The number of done and todo characters
    done=$(bc <<< "scale=0; $bar_size * $percent / 100" )
    todo=$(bc <<< "scale=0; $bar_size - $done" )

    # build the done and todo sub-bars
    done_sub_bar=$(printf "%${done}s" | tr " " "${bar_char_done}")
    todo_sub_bar=$(printf "%${todo}s" | tr " " "${bar_char_todo}")

    # output the bar
    echo -ne "\rProgress : [${done_sub_bar}${todo_sub_bar}] ${percent}%"

    if [ $total -eq $current ]; then
        echo -e "\nDONE"
    fi
}

############################
###     Index Function   ###
############################

function get_index() {
for i in "${!smiles_tab[@]}"; do
   if [[ "${smiles_tab[$i]}" = "${1}" ]]; then
       echo "${i}";
   fi
done
}

######################
### Work Directory ###
######################

work_dir="${PWD}"

DirOutput="${work_dir}/Results_Prediction/"
mkdir $DirOutput

tmp="${work_dir}/tmp/"
mkdir $tmp

#For download version of biotransformer3
#DirBiotrans="${work_dir}/biotransformer3.0jar/"
DirMetatrans="${work_dir}/MetaTrans/"
DirMetapred="${work_dir}/Meta-Predictor/"
DirCondaEnv="${work_dir}/CondaEnv/"
DirScripts="${work_dir}/Scripts/"

Script_Smiles2Smart="${DirScripts}smilestosmart.py"
Script_Smart2Smiles="${DirScripts}smarttosmiles.py"
Script_SdftoSmi="${DirScripts}sdftosmi.py"
Script_SmitoStr="${DirScripts}smitostr.py"
Script_FormulaFromSmiles="${DirScripts}formulafromsmiles.py"
Script_massFromFormula="${DirScripts}massFromFormula.py"

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
        phase1=$(echo "$option_meta" | cut -d "," -f1)
        phase2=$(echo "$option_meta" | cut -d "," -f2)

    type=$(zenity --list --title="Biotransformer 3 model" --text="Choose the type of biotransformation to use with Biotransformer3" --column="Type" --column="Description" \
    allHuman "Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step" ecbased "Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism" cyp450 "CYP450 metabolism prediction" phaseII "Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others" hgut "Human gut microbial" superbio "Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)" envimicro "Environmental microbial" )
    
    option_biotrans=$(zenity --forms --title="Biotransformer options" --text="Directly Validate to apply default values" --add-entry="The number of steps for the prediction [default=2]" --add-entry="CYP450 prediction Mode: 1=CypReact+BioTransformer rules; 2=CyProduct only; 3=CypReact+BioTransformer rules+CyProducts [Default=3]" --separator=",")
        nstep=$(echo "$option_biotrans" | cut -d "," -f1)
        cmode=$(echo "$option_biotrans" | cut -d "," -f2)

    zenity --question --title="Meta-Predictor Activation" --text="Do you want to activate Meta-Predictor ?"
    metapred_activate_answer=$?
    if [ $metapred_activate_answer -eq 0 ]; then
	    metapred_activate=true
    else [ $metapred_activate_answer -eq 1 ]
    	metapred_activate=false
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
	nstep="2"
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
if [ -z $metapred_activate ]; then
	metapred_activate=false
fi

###Keep intermediate files
if [ -z $keep_tmp ]; then
	keep_tmp=false
fi

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

if $metapred_activate; then
    if conda info --envs | grep -q metapredictor; then 
        echo "Conda environment "metapredictor" already exists"
    else 
        conda env create --name metapredictor --file ${DirCondaEnv}metapred_environment.yml
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

### Step 1 : Prepare data
echo "
     ==================================================
    ||                                                ||
    ||   Process of MetaTrans step 1 : Prepara data   ||
    ||                                                ||
     ==================================================
"

python prepare_input_file.py -input_file "${path_input}" -output_file ${outfile}

### Step 2 : Translate
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
	python translate.py -model $MODEL_FILE -src $outfile -output $OUT_FILE -n_best $BEAM -beam_size $BEAM -verbose -min_length $MIN -max_length $MAX 
done

### Step 3 : Get predictions
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
-visualise_molecules True

conda deactivate

### Déplacement des fichiers de résultats et suppression des fichiers intermédiaires
rm -r ${outfile} ${STORE} ${images}
mv ${results} ${tmp}
cd ${tmp}

#Séparation des résultats dans des fichiers distincts pour chaques molécules
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
        #n_metabo=0
        for i in ${tab_metabo[@]}; do
            #((n_metabo+=1))
            echo "${i}" >> "${MoleculeID}_Metatrans.csv"
        done
    fi
done < "Prediction_MetaTrans.csv"

######################
### META-PREDICTOR ###
######################

if ${metapred_activate}; then

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
    bash predict-top15.sh processed_data.txt ./prediction ${path_input}

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
            #SmileParent=$(echo $a | cut -d\, -f2)
            SmileMetabo=$(echo $a | cut -d\, -f3)
            read -ra tab_metabo <<< "$SmileMetabo"
            #n_metabo=0
            for i in ${tab_metabo[@]}; do
                #((n_metabo+=1))
                echo "${i}" >> "${MoleculeID}_metapred.csv"
            done
        fi
    done < "Prediction_metapred.csv"
fi

for indice in ${!tab_molecule[@]}
do
    results_file="${DirOutput}${tab_molecule[${indice}]}_CompileResults.csv"

    # Reset table
    smiles_tab=()
    formuleBrute_tab=()
    mass_tab=()
    biotrans_tab=()
    biotrans_score_tab=()
    biotrans_pathway_tab=()
    biotrans_enzyme_tab=()
    biotrans_system_tab=()
    biotrans_precursor_formule_tab=()
    biotrans_precursor_smile_tab=()
    sygma_pathway_tab=()
    sygma_score_tab=()
    sygma_tab=()
    metapred_tab=()
    metatrans_tab=()

    #########################
    ### BIOTRANSFORMERS 3 ###
    #########################

    echo "
     ===============================================================================================================
    ||                                                                                           
    ||   Process of ${tab_molecule[${indice}]} : ${tab_smiles[${indice}]} by Biotransformers3
    ||                                                                                            
     ===============================================================================================================
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
    -ocsv "${tmp}${tab_molecule[${indice}]}_Biotransformer3.csv"

    #Changement csv format
    csvformat -D ";" "${tmp}${tab_molecule[${indice}]}_Biotransformer3.csv" | gawk -v RS='"' 'NR % 2 == 0 { gsub(/\n/, "") } { printf("%s%s", $0, RT) }' > "${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut.csv"
    
    #Keep informative column
    pat_bioTrans="^InChI;InChIKey*"
    while read line
    do
        if [[ $line =~ $pat_bioTrans ]]; then
            :
        else
            info=$(echo $line | cut -d\; -f3,6,8,14,16,17,19)
            echo "${info}" >> ${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut2.csv
        fi
    done < "${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut.csv"

    sed 's/;/,/g' ${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut2.csv > ${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv

    conda activate rdkit

    while read line
    do
        smiles=$(echo "$line" | cut -d',' -f1)
        smart=$(python3 $Script_Smiles2Smart $smiles)
        new_smiles=$(python3 $Script_Smart2Smiles $smart)
        formulebrute=$(python3 $Script_FormulaFromSmiles $new_smiles)
        #formulebrute=$(echo "$line" | cut -d',' -f2)
        score=$(echo "$line" | cut -d',' -f3)
        pathway=$(echo "$line" | cut -d',' -f4)
        enzyme=$(echo "$line" | cut -d',' -f5)
        system=$(echo "$line" | cut -d',' -f6)
        smiles_precursor=$(echo "$line" | cut -d',' -f7)
        smart_precursor=$(python3 $Script_Smiles2Smart $smiles_precursor)
        new_smiles_precursor=$(python3 $Script_Smart2Smiles $smart_precursor)

        #if smiles already presents
        if [[ $(echo ${smiles_tab[@]} | fgrep -w $new_smiles) ]]; then
            index=$(get_index $new_smiles)
            biotrans_tab["${index}"]="+"
            biotrans_score_tab["${index}"]=${score}
            biotrans_pathway_tab["${index}"]=${pathway}
            biotrans_enzyme_tab["${index}"]=${enzyme}
            biotrans_system_tab["${index}"]=${system}
            biotrans_precursor_smile_tab["${index}"]=${new_smiles_precursor}
            biotrans_precursor_formule_tab["${index}"]=$(python3 $Script_FormulaFromSmiles $new_smiles_precursor)

        #if smiles not present
        else
            smiles_tab[${#smiles_tab[@]}]=${new_smiles}
            formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
            mass_tab[${#mass_tab[@]}]=$(python3 $Script_massFromFormula $formulebrute)
            index=$(get_index $new_smiles)
            biotrans_tab["${index}"]="+"
            biotrans_score_tab["${index}"]=${score}
            biotrans_pathway_tab["${index}"]=${pathway}
            biotrans_enzyme_tab["${index}"]=${enzyme}
            biotrans_system_tab["${index}"]=${system}
            biotrans_precursor_smile_tab["${index}"]=${new_smiles_precursor}
            biotrans_precursor_formule_tab["${index}"]=$(python3 $Script_FormulaFromSmiles $new_smiles_precursor)
        fi
    done < "${tmp}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv"

    ##################
    ###    SygMa   ###
    ##################

    echo "
     =================================================================================
    ||                                         
    ||   Process of ${tab_molecule[${indice}]} : ${tab_smiles[${indice}]} by Sygma
    ||                                             
     =================================================================================
    "

    pat_pathway="^>  <Pathway>"
    pat_score="^>  <Score>"

    singularity run docker://3dechem/sygma ${tab_smiles[${indice}]} -1 $phase1 -2 $phase2 >> "${tmp}${tab_molecule[${indice}]}_Sygma.sdf"

    ###Converting the sdf into smiles
    python3 ${Script_SdftoSmi} ${tmp}${tab_molecule[${indice}]}_Sygma.sdf
    mv ${PWD}/smiles.txt ${tmp}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt

    #Adding scores & pathways from the sdf file to the txt file
    while read line
    do
        if [[ $line_pre == "pathway" ]]; then
            echo $line | sed 's/,/;/g' >> "${tmp}Prediction_sygma_${tab_molecule[${indice}]}_Path.txt"
            line_pre=""
        fi

        if [[ $line_pre == "score" ]]; then
            echo $line >> "${tmp}Prediction_sygma_${tab_molecule[${indice}]}_Score.txt"
            line_pre=""
        fi

        if [[ $line =~ $pat_pathway ]]; then
            line_pre=pathway
        fi

        if [[ $line =~ $pat_score ]]; then
            line_pre=score
        fi
    done < "${tmp}${tab_molecule[${indice}]}_Sygma.sdf"

    paste -d ',' ${tmp}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt ${tmp}Prediction_sygma_${tab_molecule[${indice}]}_Path.txt ${tmp}Prediction_sygma_${tab_molecule[${indice}]}_Score.txt > ${tmp}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt

    pat_parent="^Molecule1"

    while read -r line
    do
        if [[ $line =~ $pat_parent ]]; then
            :
        else
            smiles=$(echo "$line" | cut -d',' -f2)
            pathway=$(echo "$line" | cut -d',' -f3)
            score=$(echo "$line" | cut -d',' -f4)
            smart=$(python3 $Script_Smiles2Smart $smiles)
            new_smiles=$(python3 $Script_Smart2Smiles $smart)
            
            #if smiles already presents
            if [[ $(echo ${smiles_tab[@]} | fgrep -w $new_smiles) ]]; then

                index=$(get_index $new_smiles)
                sygma_pathway_tab["${index}"]=${pathway}
                sygma_score_tab["${index}"]=${score}
                sygma_tab["${index}"]="+"
            
            #if smiles not present
            else
                smiles_tab[${#smiles_tab[@]}]=${new_smiles}
                formulebrute=$(python3 $Script_FormulaFromSmiles $new_smiles)
                formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
                mass_tab[${#mass_tab[@]}]=$(python3 $Script_massFromFormula $formulebrute)
                index=$(get_index $new_smiles)
                sygma_pathway_tab["${index}"]=${pathway}
                sygma_score_tab["${index}"]=${score}
                sygma_tab["${index}"]="+"
            fi
        fi
    done < "${tmp}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt"

    #################
    ### MetaTrans ###
    #################

    while read -r line
    do
        smiles=$(echo "$line" | cut -d',' -f1)
        smart=$(python3 $Script_Smiles2Smart $smiles)
        new_smiles=$(python3 $Script_Smart2Smiles $smart)

        #if smiles already presents
        if [[ $(echo ${smiles_tab[@]} | fgrep -w $new_smiles) ]]; then
            index=$(get_index $new_smiles)
            metatrans_tab[${index}]="+"

        #if smiles not present
        else
            smiles_tab[${#smiles_tab[@]}]=${new_smiles}
            formuleBrute_tab[${#formuleBrute_tab[@]}]=$(python3 $Script_FormulaFromSmiles $new_smiles)
            mass_tab[${#mass_tab[@]}]=$(python3 $Script_massFromFormula $formulebrute)
            index=$(get_index $new_smiles)
            metatrans_tab[${index}]="+"
        fi
    done < "${tmp}${tab_molecule[${indice}]}_Metatrans.csv"

    ######################
    ### Meta-Predictor ###
    ######################
    if ${metapred_activate}; then
        while read -r line
        do
            smiles=$(echo "$line" | cut -d',' -f2)
            smart=$(python3 $Script_Smiles2Smart $smiles)
            new_smiles=$(python3 $Script_Smart2Smiles $smart)

            #if smiles already presents
            if [[ $(echo ${smiles_tab[@]} | fgrep -w $new_smiles) ]]; then
                index=$(get_index $new_smiles)
                metapred_tab[${index}]="+"

            #if smiles not present
            else
                smiles_tab[${#smiles_tab[@]}]=${new_smiles}
                formuleBrute_tab[${#formuleBrute_tab[@]}]=$(python3 $Script_FormulaFromSmiles $new_smiles)
                mass_tab[${#mass_tab[@]}]=$(python3 $Script_massFromFormula $formulebrute)
                index=$(get_index $new_smiles)
                metapred_tab[${index}]="+"
            fi
        done < "${tmp}${tab_molecule[${indice}]}_metapred.csv"
    fi

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
    # Entete
    echo -e "FormuleBrute\tMasse(+H)\tSmiles\tSygma\tBioTransformers3\tMetaTrans\tMetaPredictor\tSygma_pathway\tBioTrans_pathway\tSygma_score\tBioTrans_score\tBioTrans_precursor\tBioTrans_precursor\tBioTrans_enzyme\tBioTrans_system\tFigure" > ${results_file}

    #Progress tool
    tasks_in_total=$( echo ${#smiles_tab[@]} )
    current_task=0

    nbmolecule=0

    for indice2 in ${!smiles_tab[@]}
    do
        ((current_task+=1))
        show_progress $current_task $tasks_in_total

        ((nbmolecule+=1))

        smiles="${smiles_tab[${indice2}]}"
        formulebrute="${formuleBrute_tab[${indice2}]}"
        masse="${mass_tab[${indice2}]}"

        sygma="${sygma_tab[${indice2}]}"
        score="${sygma_score_tab[${indice2}]}"
        pathway="${sygma_pathway_tab[${indice2}]}"

        metapred="${metapred_tab[${indice2}]}"
        metatrans="${metatrans_tab[${indice2}]}"
        
        biotrans="${biotrans_tab[${indice2}]}"
        biotrans_pathway="${biotrans_pathway_tab[${indice2}]}"
        biotrans_enzyme="${biotrans_enzyme_tab[${indice2}]}"
        biotrans_system="${biotrans_system_tab[${indice2}]}"
        biotrans_score="${biotrans_score_tab[${indice2}]}"
        biotrans_precur_for="${biotrans_precursor_formule_tab[${indice2}]}"
        biotrans_precur_smiles="${biotrans_precursor_smile_tab[${indice2}]}"

        figure="Figure_${nbmolecule}"

        echo -e "${formulebrute}\t${masse}\t${smiles}\t${sygma}\t${biotrans}\t${metatrans}\t${metapred}\t${pathway}\t${biotrans_pathway}\t${score}\t${biotrans_score}\t${biotrans_precur_for}\t${biotrans_precur_smiles}\t${biotrans_enzyme}\t${biotrans_system}\t${figure}" >> ${results_file}
        
        echo -e "Molecule${nbmolecule},${smiles}" >> "${tmp}ListeSmile.txt"
    done

    ###Structure construction
    echo "     
     =================================
    ||                               ||
    ||     Structures realisation    ||
    ||                               ||
     =================================
    "
    python3 $Script_SmitoStr -i "${tmp}ListeSmile.txt"
    mkdir "${DirOutput}${tab_molecule[${indice}]}_figures"
    mv Molecule*.jpeg "${DirOutput}${tab_molecule[${indice}]}_figures"
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