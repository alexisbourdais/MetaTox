#!/bin/bash

#======> regler probeleme sudo <======#

#####################
### Fonction HELP ###
#####################
help_msg() { 
    printf """

Usage: $0 [-h] [-i input_file] [-f format]

#############################
### Paramètre obligatoire ###
#############################

-i, --input         dossier contenant les raw

#############################
### Paramètres optionnels ###
#############################

-h, --help          show this help message and exit
-f, --format      	Format [defaut=mzXML]
"""
}

msconvert() {
######################################
### Gestion des paramètres/options ###
######################################
while [ $# -gt 0 ]; do
    key="$1"
    case $key in
        -h|--help)
            help_msg
            exit 0
            ;;
        -i|--input)
            pathBrut="$2"
            shift
            shift
            ;;
        -f|--format)
            format="$2"
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
if [ -z $pathBrut ] || [ ! -d $pathBrut ]; then
    echo "

    Erreur : Aucun input sélectionné ou le dossier $pathBrut n'existe pas !
    "
    help_msg
    exit 0
fi

###format par defaut
if [ -z $format ]; then
    format="mzXML"
fi

###dossier de sortie
pathOut="${PWD}/Fichier_mzxml/"
if [ ! -d "${pathOut}" ]; then
mkdir ${pathOut}
fi


######################
###    Program     ###
######################

echo "
###########################################################
##### Execution du programme msConvert de ProteoWizard ####
###########################################################

Les fichiers raw du dossier : ${pathBrut} 
vont être convertis au format : ${format} 
"

cd ${pathBrut}
for fichier in *.raw
do

echo "

###### Traitement du fichier ${fichier} #####

"
sudo docker run -it --rm -e WINEDEBUG=-all -v ${pathBrut}:/data \
chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert ${fichier} \
-o ${format} \
--${format} \
--inten64 \
--32 \
--filter "peakPicking vendor" \
--filter "msLevel 1-2" \
--filter "polarity positive"
done

sudo mv ${format}/*.${format} ${pathOut}
sudo rm -r ${format}/

echo "

Les fichiers au format mzXML sont dans le dossier ${pathOut}

L'execution du programme est terminée !

"
}


#filter
#Filters are applied sequentially in the order that you list them, 
#and the sequence order can make a large difference in your output. 
#In particular, the peakPicking filter must be first in line if you wish to use the vendor-supplied 
#centroiding algorithms since these use the vendor DLLs, which only operate on raw untransformed data.

#msLevel <mslevels>
#This filter selects only spectra with the indicated <mslevels>, expressed as an int_set.
