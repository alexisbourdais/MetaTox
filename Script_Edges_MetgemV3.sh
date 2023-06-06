#!/bin/bash

source ${PWD}/Tools/ProgressTool.sh

#####################
### Fonction HELP ###
#####################

help_msg() { 
    printf """

Usage: $0 [-h] [-i input_file]

#############################
### Paramètre obligatoire ###
#############################

-i, --input         fichier .csv (séparé par virgule ',')

#############################
### Paramètres optionnels ###
#############################

-h, --help          show this help message and exit

#########################
### Fichier reference ###
#########################

-> 1ère ligne   : entête (Biotransformation,DeltaMZ)
-> 1ère colonne : Biotransformation
-> 2eme colonne : DeltaMZ (0.00)
-> Séparateur   : virgule ','

Exemple:

Biotransformation,DeltaMZ
Oxidation+Demethylation+Dehydrogenation,0.03

"""
}

ID_metgem() {

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

############################
### Barre de progression ###
############################
tasks_in_total=$( nl $input | wc -l )
current_task=0

##############
### OUTPUT ###
##############
NameOutput=`basename $input .csv`_ID.csv
Directory="$(dirname "${input}")"
Output="${Directory}/${NameOutput}"

############################
### Fichier de référence ###
############################
ficcsv=${PWD}/Library/ListeBiotransformation.csv

##########################################################################
### Traitement du fichier csv + Identification des Bio-transformations ###
##########################################################################

###Pattern premiere ligne
pat='^Source*'
nbbiotrans=0

echo "
Traitement du fichier $input
"

while read a
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total
    #si premiere ligne, on copie les entetes et ajoute une colonne Shared name
    if [[ $a =~ $pat ]]; then
        echo "$(echo ${a} | cut -d\, -f1-2),DeltaMZ,Shared name,Biotransformation" > $Output
    #sinon, on copie la ligne et ajoute dans la colonne shared name source interact with target
    else
        source=$(($(echo $a | cut -d\, -f1)-1))
        target=$(($(echo $a | cut -d\, -f2)-1))
        deltamz_noRound=$(echo $a | cut -d\, -f3 | sed "s/\-//g" | awk '{printf("%.3f\n",$1)}')
        
        if (( $(echo "${deltamz_noRound} < 1" |bc -l) )); then
            deltamz=$(echo ${deltamz_noRound} | cut -c 1-4)
        #Si deltamz > 1
        else
            deltamz=$(bc -l <<< "scale = 2; ${deltamz_noRound} / 1")
        fi
    
        biotrans=$(grep -w -- "$deltamz" $ficcsv | cut -d\, -f1)
    
        if [ -n "$biotrans" ]; then
            echo "${source},${target},$(echo $a | cut -d\, -f3 | awk '{printf("%.4f\n",$1)}'),${source} (interacts with) ${target},${biotrans}" >> $Output
            ((nbbiotrans+=1))
        else
        echo "${source},${target},$(echo $a | cut -d\, -f3 | awk '{printf("%.4f\n",$1)}'),${source} (interacts with) ${target}" >> $Output
        fi
    fi
done < $input

echo "
Nombre de biotransformations détectées : $nbbiotrans

Enregistrement des résultats dans le fichier $Output

Fin de l'excecution du programme !
"
}