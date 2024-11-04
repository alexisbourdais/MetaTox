#!  /usr/bin/bash

source ${PWD}/Scripts/biotransformer3.sh
source ${PWD}/Scripts/metapredictor.sh
source ${PWD}/Scripts/metatrans.sh
source ${PWD}/Scripts/sygma.sh

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

        -m|--metapred   To activate metapredictor

        -t|--type   Type of biotransformation to use with BioTransformer3:
                       [allHuman]  : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                        ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                        cyp450     : CYP450 metabolism prediction 
                        phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                        hgut       : Human gut microbial
                        superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                        envimicro  : Environmental microbial

        -n|--nstep  The number of steps for the prediction by BioTransformers [default=1]

        -c|--cmode  CYP450 prediction Mode uses by BioTransformers: 
                        1  = CypReact+BioTransformer rules
                        2  = CyProduct only
                       [3] = CypReact+BioTransformer rules+CyProducts
                    
        -1|--phase1 Number of reaction cycles Phase 1 by SygMa [defaut=1]
        -2|--phase2 Number of reaction cycles Phase 2 by SygMa [defaut=1]

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
for i in "${!smile_tab[@]}"; do
   if [[ "${smile_tab[$i]}" = "${1}" ]]; then
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

DirScripts="${work_dir}/Scripts/"
DirBiotrans="${work_dir}/biotransformer3.0jar/"
DirMetatrans="${work_dir}/MetaTrans/"
DirMetapred="${work_dir}/Meta-Predictor/"
DirCondaEnv="${work_dir}/CondaEnv/"

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
    
    option_biotrans=$(zenity --forms --title="Biotransformer options" --text="Directly Validate to apply default values" --add-entry="The number of steps for the prediction [default=1]" --add-entry="CYP450 prediction Mode: 1=CypReact+BioTransformer rules; 2=CyProduct only; 3=CypReact+BioTransformer rules+CyProducts [Default=3]" --separator=",")
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
if [ -z $metapred_activate ]; then
	metapred_activate=false
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

###############
### Program ###
###############

for indice in ${!tab_molecule[@]}
do
    sygma_function --script ${DirScripts} --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput} -1 ${phase1} -2 ${phase2}
    biotransformer_function --script ${DirScripts} --molecule ${tab_molecule[${indice}]} --smile ${tab_smiles[${indice}]} --outdir ${DirOutput} --biodir ${DirBiotrans} --type ${type} --nstep ${nstep} --cmode ${cmode}
done

metatrans_function --script ${DirScripts} --input "${path_input}" --outdir ${DirOutput} --metadir ${DirMetatrans} --min ${MIN} --max ${MAX} --beam ${BEAM}

if ${metapred_activate}; then
    metapredictor_function --script ${DirScripts} --input "${path_input}" --outdir ${DirOutput} --metadir ${DirMetapred}
fi

#######################
### Compile results ###
#######################

for indice in ${!tab_molecule[@]}
do
    molecule="${tab_molecule[${indice}]}"

    sygma_file="$DirOutput/${molecule}_Sygma.csv"
    trans_file="$DirOutput/${molecule}_Metatrans.csv"
    pred_file="$DirOutput/${molecule}_metapred.csv"

    pattern="^Métabolites*"
    results_file="${tab_molecule[${indice}]}_CompileResults.csv"

    # Table
    smile_tab=()
    sygma_pathway_tab=()
    sygma_score_tab=()
    formuleBrute_tab=()
    mass_tab=()
    metapred_tab=()
    metapred_fig_tab=()
    metatrans_tab=()
    metatrans_fig_tab=()
    sygma_tab=()
    sygma_fig_tab=()

    #############
    ### Sygma ###
    #############

    echo "
    Compile $molecule results"

    tasks_sygma=$( nl "$sygma_file" | wc -l )
    tasks_pred=$( nl "$pred_file" | wc -l )
    tasks_trans=$( nl "$trans_file" | wc -l )

    tasks_in_total=$((tasks_sygma+tasks_pred+tasks_trans))
    current_task=0

    while read -r line
    do
        ((current_task+=1))
        show_progress $current_task $tasks_in_total

        if [[ $line =~ $pattern ]]; then
            :

        else
            id=$(echo "$line" | cut -d',' -f1 | sed -e 's/Molecule//')
            smile=$(echo "$line" | cut -d',' -f2)
            pathway=$(echo "$line" | cut -d',' -f3)
            score=$(echo "$line" | cut -d',' -f4)
            formulebrute=$(echo "$line" | cut -d',' -f5)
            masse=$(echo "$line" | cut -d',' -f6)

            smile_tab[${#smile_tab[@]}]=${smile}
            sygma_pathway_tab[${#sygma_pathway_tab[@]}]=${pathway}
            sygma_score_tab[${#sygma_score_tab[@]}]=${score}
            formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
            mass_tab[${#mass_tab[@]}]=${masse}
            index=$(get_index $smile)
            sygma_tab[${index}]="+"
            sygma_fig_tab[${index}]=${id}
        fi
    done < "$sygma_file"
    rm $sygma_file

    #################
    ### MetaTrans ###
    #################

    while read -r line
    do
        ((current_task+=1))
        show_progress $current_task $tasks_in_total

        if [[ $line =~ $pattern ]]; then
            :

        else
            id=$(echo "$line" | cut -d',' -f1 | sed -e 's/Metabolite//')
            smile=$(echo "$line" | cut -d',' -f2)
            formulebrute=$(echo "$line" | cut -d',' -f3)
            masse=$(echo "$line" | cut -d',' -f4)

            #if smile already presents
            if [[ $(echo ${smile_tab[@]} | fgrep -w $smile) ]]; then
                index=$(get_index $smile)
                metatrans_tab[${index}]="+"
                metatrans_fig_tab[${index}]=${id}

            #if smile not present
            else
                smile_tab[${#smile_tab[@]}]=${smile}
                formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
                mass_tab[${#mass_tab[@]}]=${masse}
                index=$(get_index $smile)
                metatrans_tab[${index}]="+"
                metatrans_fig_tab[${index}]=${id}
            fi
        fi

    done < "$trans_file"
    rm $trans_file

    ######################
    ### Meta-Predictor ###
    ######################

    while read -r line
    do
        ((current_task+=1))
        show_progress $current_task $tasks_in_total

        if [[ $line =~ $pattern ]]; then
            :

        else
            id=$(echo "$line" | cut -d',' -f1 | sed -e 's/Metabolite//')
            smile=$(echo "$line" | cut -d',' -f2)
            formulebrute=$(echo "$line" | cut -d',' -f3)
            masse=$(echo "$line" | cut -d',' -f4)

            #if smile already presents
            if [[ $(echo ${smile_tab[@]} | fgrep -w $smile) ]]; then
                index=$(get_index $smile)
                metapred_tab[${index}]="+"
                metapred_fig_tab[${index}]=${id}

            #if smile not present
            else
                smile_tab[${#smile_tab[@]}]=${smile}
                formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
                mass_tab[${#mass_tab[@]}]=${masse}
                index=$(get_index $smile)
                metapred_tab[${index}]="+"
                metapred_fig_tab[${index}]=${id}
            fi
        fi
    done < "$pred_file"
    rm $pred_file

    ###################
    ### Comparaison ###
    ###################

    # Entete
    echo -e "FormuleBrute\tMasse(+H)\tSmile\tSygma_pathway\tSygma_score\tSygma\tMetaTrans\tMetaPredictor\tSygMa-Figure\tMetaTrans-Figure\tMetaPred-Figure" > ${results_file}

    for indice in ${!smile_tab[@]}
    do
        smile="${smile_tab[${indice}]}"
        formulebrute="${formuleBrute_tab[${indice}]}"
        masse="${mass_tab[${indice}]}"
        score="${sygma_score_tab[${indice}]}"
        pathway="${sygma_pathway_tab[${indice}]}"

        metapred="${metapred_tab[${indice}]}"
        metapred_fig="${metapred_fig_tab[${indice}]}"

        metatrans="${metatrans_tab[${indice}]}"
        metatrans_fig="${metatrans_fig_tab[${indice}]}"

        sygma="${sygma_tab[${indice}]}"
        sygma_fig="${sygma_fig_tab[${indice}]}"

        echo -e "${formulebrute}\t${masse}\t${smile}\t${pathway}\t${score}\t${sygma}\t${metatrans}\t${metapred}\t${sygma_fig}\t${metatrans_fig}\t${metapred_fig}" >> ${results_file}
    done
done

echo "
Recording results in : ${DirOutput}

Execution completed !
"