#!  /usr/bin/bash

metatrans_function () {

#############
### Input ###
#############

while [ $# -gt 0 ] ; do
    key="$1"
    case $key in
        -h|--help)
            help_msg
            exit 0
            ;;
        -i|--input)
            input="$2"
            shift
            shift
            ;;
        -o|--outdir)
            DirOutput="$2"
            shift
            shift
            ;;
        -d|--metadir)
            DirMetatrans="$2"
            shift
            shift
            ;;
        -m|--min)
            MIN="$2"
            shift
            shift
            ;;
        -M|--max)
            MAX="$2"
            shift
            shift
            ;;
        -b|--beam)
            BEAM="$2"
            shift
            shift
            ;;
        -r|--script)
            DirScripts="$2"
            shift
            shift
            ;;
        *)
            echo "Argument non défini : '$key'"
            exit 1                              
    esac
done

#############
### Utils ###
#############

Script_SdftoSmi="${DirScripts}sdftosmi.py"
Script_SmitoStr="${DirScripts}smitostr.py"
Script_FormulaFromSmiles="${DirScripts}formulafromsmiles.py"
Script_massFromFormula="${DirScripts}massFromFormula.py"

###############
### Program ###
###############

DirOutputFig=${DirOutput}MetaTrans_Figures
mkdir $DirOutputFig

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

python prepare_input_file.py -input_file ${input} -output_file ${outfile}

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
-input_file ${input} \
-output_file ${results} \
-beam_size ${BEAM} \
-visualise_molecules True

conda deactivate

### Déplacement des fichiers de résultats et suppression des fichiers intermédiaires
rm ${outfile}
rm -r ${STORE}
mv ${results} ${DirOutput}
mv ${images} ${DirOutputFig}

### Retraitement fichier de résultats MetaTrans ###
cd ${DirOutput}
pat='^Molecule'

#Séparation des résultats dans des fichiers distincts pour chaques molécules
while read a
do
    if [[ $a =~ $pat ]]; then
        :

    else
        MoleculeID=$(echo $a | cut -d\, -f1)
        SmileParent=$(echo $a | cut -d\, -f2)
        SmileMetabo=$(echo $a | cut -d\, -f3)
        read -ra tab_metabo <<< "$SmileMetabo"
        n_metabo=0
        for i in ${tab_metabo[@]}; do
            ((n_metabo+=1))
            echo "Metabolite${n_metabo},${i}" >> "${MoleculeID}_MetatransV1.txt"
        done
    fi
done < "Prediction_MetaTrans.csv"

rm ${DirOutput}${results}

conda activate rdkit

for indice in ${!tab_molecule[@]}
do
    #Calcul Mass+H from brut formula
    while read line
    do
        smiles=$(echo $line | cut -d\, -f2)
        python3 $Script_FormulaFromSmiles $smiles
    done < ${tab_molecule[${indice}]}_MetatransV1.txt >> ${tab_molecule[${indice}]}_Metatrans_BrutFormula.txt

    #Calcul Mass+H from brut formula
    while read line
    do
        formula=$(echo $line)
        python3 $Script_massFromFormula $formula
    done < ${tab_molecule[${indice}]}_Metatrans_BrutFormula.txt >> ${tab_molecule[${indice}]}_Metatrans_Mass.txt

    #Complilation des colonnes
    paste -d ',' ${tab_molecule[${indice}]}_MetatransV1.txt ${tab_molecule[${indice}]}_Metatrans_BrutFormula.txt ${tab_molecule[${indice}]}_Metatrans_Mass.txt > ${tab_molecule[${indice}]}_MetatransVF.txt

    #Création des entêtes
    echo "Métabolites,SMILES,FormuleBrute,Masse(+H)" > ${tab_molecule[${indice}]}_Metatrans.csv
    while read line
    do
        echo $line
    done < ${tab_molecule[${indice}]}_MetatransVF.txt >> ${tab_molecule[${indice}]}_Metatrans.csv

    rm ${tab_molecule[${indice}]}_MetatransV1.txt ${tab_molecule[${indice}]}_Metatrans_BrutFormula.txt ${tab_molecule[${indice}]}_Metatrans_Mass.txt ${tab_molecule[${indice}]}_MetatransVF.txt
done

conda deactivate
}