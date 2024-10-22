#!  /usr/bin/bash

metapredictor_function () {

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
            DirMetapred="$2"
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


DirOutputFig=${DirOutput}MetaPred_Figures
mkdir $DirOutputFig

conda activate metapredictor

cd ${DirMetapred}

echo "
    #####   Process of MetaPredictor   #####
"

python prepare_input_file.py -input_file $input -output_file processed_data.txt
./predict-top15.sh  processed_data.txt  ./prediction  $input

conda deactivate

### Déplacement des fichiers de résultats et suppression des fichiers intermédiaires
mv prediction/predict.csv ${DirOutput}Prediction_metapred.csv
mv Figures/* ${DirOutputFig}
rm prediction/*
rm processed_data.txt

### Retraitement fichier de résultats metapred ###
cd ${DirOutput}
pat='^Name'

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
            echo "Metabolite${n_metabo},${i}" >> "${MoleculeID}_metapredV1.txt"
        done
    fi
done < "Prediction_metapred.csv"

rm Prediction_metapred.csv

conda activate my-rdkit-env

for indice in ${!tab_molecule[@]}
do
    #Calcul Mass+H from brut formula
    while read line
    do
        smiles=$(echo $line | cut -d\, -f2)
        python3 $Script_FormulaFromSmiles $smiles
    done < ${tab_molecule[${indice}]}_metapredV1.txt >> ${tab_molecule[${indice}]}_metapred_BrutFormula.txt

    #Calcul Mass+H from brut formula
    while read line
    do
        formula=$(echo $line)
        python3 $Script_massFromFormula $formula
    done < ${tab_molecule[${indice}]}_metapred_BrutFormula.txt >> ${tab_molecule[${indice}]}_metapred_Mass.txt

    #Complilation des colonnes
    paste -d ',' ${tab_molecule[${indice}]}_metapredV1.txt ${tab_molecule[${indice}]}_metapred_BrutFormula.txt ${tab_molecule[${indice}]}_metapred_Mass.txt > ${tab_molecule[${indice}]}_metapredVF.txt

    #Création des entêtes
    echo "Métabolites,SMILES,FormuleBrute,Masse(+H)" > ${tab_molecule[${indice}]}_metapred.csv
    while read line
    do
        echo $line
    done < ${tab_molecule[${indice}]}_metapredVF.txt >> ${tab_molecule[${indice}]}_metapred.csv

    rm ${tab_molecule[${indice}]}_metapredV1.txt ${tab_molecule[${indice}]}_metapred_BrutFormula.txt ${tab_molecule[${indice}]}_metapred_Mass.txt ${tab_molecule[${indice}]}_metapredVF.txt
done

conda deactivate
}