#!  /usr/bin/bash

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

-i, --input         fichier .selfloop

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


ID_biotransformation() {

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

################
### DOS2UNIX ###
################
REQUIRED_PKG="dos2unix"
PKG_OK=$(dpkg-query -W --showformat='${Status}\n' $REQUIRED_PKG|grep "install ok installed")
echo Checking for $REQUIRED_PKG: $PKG_OK
if [ "" = "$PKG_OK" ]; then
  echo "No $REQUIRED_PKG. Setting up $REQUIRED_PKG."
  sudo apt-get --yes install $REQUIRED_PKG
fi

file $input | grep CRLF && dos2unix $input


############################
### Barre de progression ###
############################
tasks_in_total=$( nl $input | wc -l )
current_task=0

##############
### OUTPUT ###
##############
NameOutput=`basename $input .selfloop`_Edges_ID.selfloop
Directory="$(dirname "${input}")"
Output="${Directory}/${NameOutput}"

############################
### Fichier de référence ###
############################
ficcsv=${PWD}/Library/ListeBiotransformation.csv


##############################################
### Identification des Bio-transformations ###
##############################################
echo "
Identification des Bio-transformations
"
nbbiotrans=0

###Pattern premiere ligne fichier
pat='^CLUSTERID1'

while read a
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total
    #si premiere ligne, on copie les entetes et ajoute une colonne Biotransformation
    if [[ $a =~ $pat ]]; then
        printf '%s\t%s\n' "${a}" "Biotransformation" > $Output

    else
        #deltamz = variable de la colonne deltamz arrondie au 3 chiffre, sans le moins(-)
        deltamz_noRound=$(echo $a | cut -d' ' -f3 | sed "s/\-//g" | awk '{printf("%.3f\n",$1)}')
        #Si deltamz < 1
        if (( $(echo "${deltamz_noRound} < 1" |bc -l) )); then
            deltamz=$(echo ${deltamz_noRound} | cut -c 1-4)
        #Si deltamz > 1
        else
            deltamz=$(bc -l <<< "scale = 2; ${deltamz_noRound} / 1")
        fi
        #recherche de la biotransformation correspondante dans le fichier reference(champs separe par ',')
        biotrans=$(grep -w -- "$deltamz" $ficcsv | cut -d\, -f1)
        #annotation des liens fait par GNPS
        edgeAnnot=$(echo $a | cut -d' ' -f8)

        if [ -n "$biotrans" -a -z "$edgeAnnot" ]; then
            ((nbbiotrans+=1))
            echo -e "${a}\t\t${biotrans}" >> $Output

        #si biotransformation correspondante, copie la ligne et ajoute la biotrans dans la dernire colonne
        elif [ -n "$biotrans" -a -n "$edgeAnnot" ]; then
            ((nbbiotrans+=1))
            echo -e "${a}\t${biotrans}" >> $Output

        #si juste annotation de GNPS, copie la ligne
        else
            echo -e "$a" >> $Output
        fi
    fi
done < $input

echo "
Nombre de biotransformations détectées : $nbbiotrans

Enregistrement des résultats dans le fichier ${Output}

Fin de l'excecution du programme !
"
}