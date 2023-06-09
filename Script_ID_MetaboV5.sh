#!  /usr/bin/bash

#===> ajoutez fin programme si erreur dans fichier de reference<===#

source ${PWD}/Tools/ProgressTool.sh

#####################
### Fonction HELP ###
#####################
help_msg() { 
    printf """

Usage: $0 [-h] [-i input_file] [-m float] [-t float] [-r ref_file] [-p string]

#############################
### Paramètre obligatoire ###
#############################

-i, --input         fichier .clustersummary

#############################
### Paramètres optionnels ###
#############################

-h, --help          show this help message and exit
-m, --mass_pre      precision masse [defaut = 0.001]
-t, --time_pre      precision RT [defaut = 0.1]
-r, --ref           fichier csv de reference [defaut = ListeMetabolitesEutylone.csv]
-p, --pat           pattern 1ere ligne du fichier de reference [defaut = Metabolites]

#########################
### Fichier reference ###
#########################

-> 1ère ligne   : entête (Metabolites,Masse,RT)
-> 1ère colonne : nom du métabolites
-> 2eme colonne : masse (0.000)
-> 3ème colonne : RT (0.00)
-> Séparateur   : virgule ','

Exemple:

Metabolites,Masse,RT
MoleculeName,0.000,0.00

"""
}

ID_metabo() {

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
        -m|--mass_pre)
            precision_mass="$2"
            shift
            shift
            ;;
        -t|--time_pre)
            precision_rt="$2"
            shift
            shift
            ;;
        -r|--ref)
            ficcsv="$2"
            shift
            shift
            ;;
        -p|--pat)
            pat2="^$2"
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

###Si probleme reference
if [ ! -f "$ficcsv" ]; then
    ficcsv=./ListeMetabolitesEutylone.csv
    if [ ! -f "$ficcsv" ]; then
        echo "
        Erreur : Le fichier de référence $ficcsv n'existe pas !
        "
        exit 0
    fi
fi

###Précision par défaut
listefloat="^[0-9]+([.][0-9]+)?$"
if [ -z $precision_mass ] || ! [[ $precision_mass =~ $listefloat ]]; then
    precision_mass=0.001
fi

if [ -z $precision_rt ] || ! [[ $precision_rt =~ $listefloat ]]; then
    precision_rt=0.1
fi

###Pattern fichier de reference
if [ -z $pat2 ]; then
    pat2='^Metabolites'
fi

echo "
*** Paramètres sélectionnés ***

- input              : $input
- precision_mass     : $precision_mass
- precision_rt       : $precision_rt
- ref_file           : $ficcsv
- pattern ref_file   : $(echo $pat2 | sed "s/\^//g")
"

############################
### Barre de progression ###
############################
tasks_in_total=$( nl $input | wc -l )
current_task=0


##############
### OUTPUT ###
##############
NameOutput=`basename $input .clustersummary`_Metabo_ID.clustersummary
Directory="$(dirname "${input}")"
Output="${Directory}/${NameOutput}"


##############################################
### Enregistrement du fichier de référence ###
##############################################
declare -a tab_metabo
declare -a tab_mass_inf
declare -a tab_mass_sup
declare -a tab_rt_inf
declare -a tab_rt_sup

while read LigneRef
do
    if [[ $LigneRef =~ $pat2 ]]; then
        :

    else
        MetaboRef=$(echo $LigneRef | cut -d\, -f1)
        tab_metabo[${#tab_metabo[@]}]=${MetaboRef}

        MassRef=$(echo $LigneRef | cut -d\, -f2 | awk '{printf("%.3f\n",$1)}')
        MassRef_inf=$(echo "${MassRef}-${precision_mass}" |bc )
        MassRef_sup=$(echo "${MassRef}+${precision_mass}" |bc )
        tab_mass_inf[${#tab_mass_inf[@]}]=${MassRef_inf}
        tab_mass_sup[${#tab_mass_sup[@]}]=${MassRef_sup}

        RTRef=$(echo $LigneRef | cut -d\, -f3)
        RTRef_inf=$(echo "${RTRef}-${precision_rt}" |bc )
        RTRef_sup=$(echo "${RTRef}+${precision_rt}" |bc )
        tab_rt_inf[${#tab_rt_inf[@]}]=${RTRef_inf}
        tab_rt_sup[${#tab_rt_sup[@]}]=${RTRef_sup}
    fi
done < $ficcsv

echo "
Fin de l'enregistrement du fichier de référence

Identification des métabolites
"


######################################
### Identification des métabolites ###
######################################

###Pattern Première Ligne
pat='^AllGroups'
nbhit=0
while read a
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total

    if [[ $a =~ $pat ]]; then
        echo -e "${a}\t"Identification"" > $Output

    else
        Match=False
        Mass=$(echo $a | cut -d' ' -f25 | awk '{printf("%.3f\n",$1)}')
        RT=$(echo $a | cut -d' ' -f16 | awk '{printf("%.3f\n",$1)}')

        for indice in ${!tab_metabo[@]}
        do
            if [ $(echo "(${RT} >= ${tab_rt_inf[${indice}]}) && (${RT} <= ${tab_rt_sup[${indice}]}) && (${Mass} >= ${tab_mass_inf[${indice}]}) && (${Mass} <= ${tab_mass_sup[${indice}]})" |bc -l) -eq 1 ]; then 
                echo -e "\t${a}\t${tab_metabo[${indice}]}" >> $Output
                echo "      Identification de ${tab_metabo[${indice}]}; Masse = ${Mass} ; RT = ${RT}"
                ((nbhit+=1))
                Match=True
                break
            fi
        done

        if [ $Match == False ]; then
            echo -e "\t${a}" >> $Output
        fi
    fi
done < $input

echo "
Enregistrement des résultats dans le fichier ${Output}
Nombre de métabolites détectés : $nbhit

Fin de l'excecution du programme !
"
}
