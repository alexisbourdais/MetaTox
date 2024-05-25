#!/bin/bash

#===> installation miniconda, java, docker si non installés <===#


#####################
### Fonction HELP ###
#####################

help_msg() { 
    printf """

Usage: $0 [-h] [-i input_file] [-m float] [-t float] [-r ref_file] [-p string]

#############################
### Paramètre obligatoire ###
#############################

-i, --input         fichier .txt

#############################
### Paramètres optionnels ###
#############################

-h, --help      Show this help message and exit
-b, --beam      Top métabolites [defaut = 5]
-m, --min       Taille minimum des métabolites (en SMILE) [defaut = 5]
-M, --max       Taille maximum des métabolites (en SMILE) [defaut = 120]
-t, --type      Type Biotransformers3 [defaut = allhuman] (ecbased,cyp450,phaseII,hgut,superbio,envimicro)

#####################
### Fichier texte ###
#####################

-> 1ère colonne : ID/Nom molecule
-> 2eme colonne : code SMILE
-> Séparateur   : virgule ','

Exemple:

Eutylone,CCC(C(=O)C1=CC2=C(C=C1)OCO2)NCC

#################
### Prérequis ###
#################

- Necessite Docker pour l'outil Sygma (https://docs.docker.com/)
- Necessite Conda pour l'outil MetaTrans (https://docs.conda.io/projects/conda/en/stable/user-guide/index.html)
- Necessite Java pour l'outil BioTransformers3 (https://www.java.com/fr/)
"""
}

predict_metabo() {

######################################
### Gestion des paramètres/options ###
######################################

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
        -t|--type)
            type="$2"
            shift
            shift
            ;;
        -b|--beam)
            BEAM="$2"
            shift
            shift
            ;;
        *)
            echo "Argument non défini : '$key'"
            exit 1                              
    esac
done

###Si probleme d'input
if [ -z $input ] || [ ! -f $input ]; then
    echo "

    Erreur : Aucun input sélectionné ou le fichier $input n'existe pas !
    "
    help_msg
    exit 0
fi

###Mode par défaut
if [ -z $type ]; then
	type="allHuman"
fi

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


echo "
*** Paramètres actuels ***

- MetaTrans:
BEAM = ${BEAM} (5: top 10 metabolites, 10: top 20 métabolites)
MIN = ${MIN}
MAX = ${MAX}

- Biotransformers3:
Type = ${type}
"

######################
### Work Directory ###
######################

DirOutput="${PWD}/Results_Prediction/"
mkdir $DirOutput
DirOutputFig=${DirOutput}MetaTrans_Figures
mkdir $DirOutputFig
DirFigSygma=${DirOutput}Sygma_Figures
mkdir $DirFigSygma
DirFigBioTrans=${DirOutput}BioTrans_Figures
mkdir $DirFigBioTrans

DirBiotrans="${PWD}/Tools/Biotransformer3.0/"
DirMetatrans="${PWD}/Tools/MetaTrans-master/"
Script_SdftoSmi="${PWD}/Tools/sdftosmi.py"
Script_SmitoStr="${PWD}/Tools/smitostr.py"
Script_SmitoStr2="${PWD}/Tools/smitostr_2.py"
Script_FormulaFromSmiles="${PWD}/Tools/formulafromsmiles.py"
Script_massFromFormula="${PWD}/Tools/massFromFormula.py"

###############
###  Input  ###
###############

REQUIRED_PKG="dos2unix"
PKG_OK=$(dpkg-query -W --showformat='${Status}\n' $REQUIRED_PKG|grep "install ok installed")
echo Checking for $REQUIRED_PKG: $PKG_OK
if [ "" = "$PKG_OK" ]; then
  echo "No $REQUIRED_PKG. Setting up $REQUIRED_PKG."
  sudo apt-get --yes install $REQUIRED_PKG
fi

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

#########################
### Environment conda ###
#########################

eval "$(conda shell.bash hook)"
#conda init --all

if conda info --envs | grep -q metatrans; then 
echo "Conda environment "metatrans" already exists"
else 
conda env create --name metatrans --file ./Tools/metatrans_environment.yml
fi

if conda info --envs | grep -q my-rdkit-env; then 
echo "Conda environment "my-rdkit-env" already exists"
else 
conda env create --name my-rdkit-env --file ./Tools/rdkit_environment.yml
fi

##############
#### Sygma ###
##############

sudo service docker start
conda activate my-rdkit-env

pat_pathway="^>  <Pathway>"
pat_score="^>  <Score>"

for indice in ${!tab_molecule[@]}
do
    echo "
    #####   Process of ${tab_molecule[${indice}]} : ${tab_smiles[${indice}]} by Sygma   #####
    "
    sudo docker run 3dechem/sygma ${tab_smiles[${indice}]} >> "${DirOutput}${tab_molecule[${indice}]}_Sygma.sdf"

    #--phase1 Number of phase 1 cycles (default: 1)
    #--phase2 Number of phase 2 cycles (default: 1)
    #-o, --outputtype Molecule output type (default: sdf)

    ###Conversion des sdf en smiles
    python3 ${Script_SdftoSmi} ${DirOutput}${tab_molecule[${indice}]}_Sygma.sdf
    mv ${PWD}/smiles.txt ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt

    ###Réalisation des structures
    python3 $Script_SmitoStr -i ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt
    mkdir ${DirFigSygma}/${tab_molecule[${indice}]}
    mv Molecule*.jpeg ${DirFigSygma}/${tab_molecule[${indice}]}

    #Ajout des scores & pathway du fichier sdf dans le fichier txt
    while read line
    do
        if [[ $line_pre == "pathway" ]]; then
            echo $line | sed 's/,/;/g' >> "${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Path.txt"
            line_pre=""
        fi

        if [[ $line_pre == "score" ]]; then
            echo $line >> "${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Score.txt"
            line_pre=""
        fi

        if [[ $line =~ $pat_pathway ]]; then
            line_pre=pathway
        fi

        if [[ $line =~ $pat_score ]]; then
            line_pre=score
        fi
    done < "${DirOutput}${tab_molecule[${indice}]}_Sygma.sdf"

    paste -d ',' ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Path.txt ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Score.txt > ${DirOutput}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt
    rm ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_smiles.txt ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Path.txt ${DirOutput}Prediction_sygma_${tab_molecule[${indice}]}_Score.txt

    #Convert smiles into brut formula
    while read line
    do
        smiles=$(echo $line | cut -d\, -f2)
        python3 $Script_FormulaFromSmiles $smiles
    done < ${DirOutput}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt >> ${DirOutput}${tab_molecule[${indice}]}_Sygma_formuleBrute.txt

    #Calcul Mass+H from brut formula
    while read line
    do
        formula=$(echo $line)
        python3 $Script_massFromFormula $formula
    done < ${DirOutput}${tab_molecule[${indice}]}_Sygma_formuleBrute.txt >> ${DirOutput}${tab_molecule[${indice}]}_Sygma_mass.txt

    #Compilation des colonnes
    paste -d ',' ${DirOutput}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt ${DirOutput}${tab_molecule[${indice}]}_Sygma_formuleBrute.txt ${DirOutput}${tab_molecule[${indice}]}_Sygma_mass.txt > ${DirOutput}${tab_molecule[${indice}]}_SygmaVF.txt
    rm ${DirOutput}${tab_molecule[${indice}]}_Sygma_SmilesPathScore.txt ${DirOutput}${tab_molecule[${indice}]}_Sygma_formuleBrute.txt ${DirOutput}${tab_molecule[${indice}]}_Sygma_mass.txt 

    #Création des entêtes
    echo "Métabolites,SMILES,Pathway,Score,FormuleBrute,Masse(+H)" > ${DirOutput}${tab_molecule[${indice}]}_Sygma.csv
    while read line
    do 
        echo $line | sed "s/;//g"
    done < ${DirOutput}${tab_molecule[${indice}]}_SygmaVF.txt >> ${DirOutput}${tab_molecule[${indice}]}_Sygma.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_SygmaVF.txt
done

########################
### Bio-transformers ###
########################

cd ${DirBiotrans}

pat_bioTrans="^InChI;InChIKey"

for indice in ${!tab_molecule[@]}
do
    echo "
    #####   Process of ${tab_molecule[${indice}]} : ${tab_smiles[${indice}]} by Biotransformers3   #####
    "

    java -jar BioTransformer3.0.jar \
    -b ${type} \
    -k pred \
    -cm 3 \
    -ismi "${tab_smiles[${indice}]}" \
    -ocsv "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3.csv"

    #Changement de format du csv
    csvformat -D ";" "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3.csv" | gawk -v RS='"' 'NR % 2 == 0 { gsub(/\n/, "") } { printf("%s%s", $0, RT) }' > "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut.csv"

    #Suppression des colonnes inutiles
    molecule_nb=0
    while read line
    do
        if [[ $line =~ $pat_bioTrans ]]; then
            :
        else
            ((molecule_nb+=1))
            info=$(echo $line | cut -d\; -f3,6,8,14,16,17)
            echo "Molecule_${molecule_nb};${info}" >> ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut2.csv
        fi
    done < ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut.csv
    sed 's/;/,/g' ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut2.csv > ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut2.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3.csv

    ###Réalisation des structures
    python3 $Script_SmitoStr -i "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv"
    mkdir ${DirFigBioTrans}/${tab_molecule[${indice}]}
    mv Molecule*.jpeg ${DirFigBioTrans}/${tab_molecule[${indice}]}

    #Calcul Mass+H from brut formula
    while read line
    do
        formula=$(echo $line | cut -d\, -f3)
        python3 $Script_massFromFormula $formula
    done < "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv" >> "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv"

    #Fusion des colonnes
    paste -d ',' ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv > ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv

    #Creation des entetes
    echo "Métabolites,SMILES,FormuleBrute,ALogP,Pathway,Enzyme,Systeme,Masse(+H)" > ${DirOutput}${tab_molecule[${indice}]}_BioTransformer3.csv
    while read line
    do 
        echo $line
    done < "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv" >> ${DirOutput}${tab_molecule[${indice}]}_BioTransformer3.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv

done

conda deactivate

#################
### MetaTrans ###
#################

conda activate metatrans
cd $DirMetatrans

#Variables reprises de MetaTrans
outfile=processed_data.txt
results=Prediction_MetaTrans.csv
images=./Figures/*
STORE=predictions/  #directory for output process
mkdir $STORE

### Step 1 : Prepare data
echo "
#####		Process of MetaTrans step 1 : Prepara data		#####
"
python prepare_input_file.py -input_file ${input} -output_file ${outfile}

### Step 2 : Translate
echo "
#####		Process of MetaTrans step 2 : Translate		#####
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
#####		Process of MetaTrans step 2 : Get prediction		#####
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

conda activate my-rdkit-env

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

echo "
Enregistrement des résultats dans le dossier ${DirOutput}

Fin de l'excecution du programme !
"
}
