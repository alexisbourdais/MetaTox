#!/bin/bash

#===> ajoutez installation miniconda, java, docker si non installé <===#


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

-h, --help          show this help message and exit
-b, --beam      	Top métabolites [defaut = 5]
-m, --min			taille minimum des métabolites (en SMILE) [defaut = 5]
-M, --max			taille maximum des métabolites (en SMILE) [defaut = 120]
-t, --type			type Biotransformers3 [defaut = allhuman] (ecbased,cyp450,phaseII,hgut,superbio,envimicro)

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

######################
### Gestion erreur ###
######################
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

ListeSmileMolecule="${input}"

DirOutput="${PWD}/Results_Prediction/"
mkdir $DirOutput
DirOutputFig=${DirOutput}MetaTrans_Figures
mkdir $DirOutputFig

DirBiotrans="${PWD}/Tools/Biotransformer3.0/"
DirMetatrans="${PWD}/Tools/MetaTrans-master/"


#########################
### Environment conda ###
#########################

eval "$(conda shell.bash hook)"
if conda info --envs | grep -q metatrans; then 
echo "Conda environment "metatrans" already exists"
else 
conda env create --name metatrans --file environment_metatrans.yml
fi


##############
#### Sygma ###
##############

#--phase1 Number of phase 1 cycles (default: 1)
#--phase2 Number of phase 2 cycles (default: 1)
#-o, --outputtype Molecule output type (default: sdf)

while read a
do
Molecule=$(echo $a | cut -d',' -f1)
Smile=$(echo $a | cut -d',' -f2)
echo "
#####		Process of ${Molecule} : ${Smile} by Sygma		#####
"
sudo docker run 3dechem/sygma ${Smile} >> "${DirOutput}Prediction_sygma_${Molecule}.sdf"
done < $ListeSmileMolecule


########################
### Bio-transformers ###
########################

cd ${DirBiotrans}

while read a
do
Molecule=$(echo $a | cut -d',' -f1)
Smile=$(echo $a | cut -d',' -f2)
echo "
#####		Process of ${Molecule} : ${Smile} by Biotransformer3	#####
"
java -jar BioTransformer3.0.jar \
-b ${type} \
-k pred \
-cm 3 \
-ismi "${Smile}" \
-ocsv "${DirOutput}Prediction_biotransformers_${Molecule}.csv"
done < $ListeSmileMolecule


#################
### MetaTrans ###
#################

conda activate metatrans
cd $DirMetatrans

outfile=processed_data.txt
results=Prediction_MetaTrans.csv
images=./Figures/*
STORE=predictions/  #directory for output process
mkdir $STORE

### Step 1 : Prepare data
echo "
#####		Process of MetaTrans step 1 : Prepara data		#####
"
python prepare_input_file.py -input_file ${ListeSmileMolecule} -output_file ${outfile}

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
-input_file ${ListeSmileMolecule} \
-output_file ${results} \
-beam_size ${BEAM} \
-visualise_molecules True

### Déplacement des fichiers de résultats et suppression des fichiers intermédiaires
rm ${outfile}
rm -r ${STORE}
mv ${results} ${DirOutput}
mv ${images} ${DirOutputFig}

echo "
Enregistrement des résultats dans le fichier ${DirOutput}

Fin de l'excecution du programme !
"
}