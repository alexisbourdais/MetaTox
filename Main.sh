#!  /usr/bin/bash

source ${PWD}/Tools/Script_ID_MetaboV5.sh
source ${PWD}/Tools/Script_ID_EdgesV4.sh
source ${PWD}/Tools/Script_Prediction_MetaboV2.sh
source ${PWD}/Tools/Script_Edges_MetgemV3.sh
source ${PWD}/Tools/pwizV2.sh


zenity --info --text "
Lancement du programme
"

mode=$(zenity --list --title="Choix du mode" --text="Indiquez quel mode utiliser" --column="Mode" --column="Description" \
Metabolites "Identifier les métabolites d'un fichier .clustersummary à partir d'un fichier de référence" \
Biotransformations "Identifier les Bio-transformations d'un fichier .selfloop" \
InSilico "Prédiction in silico par Bio-transformer3, Sygma et MetaTrans" \
MetGem "Traitement du fichier csv contenant les informations des liens pour visualisation dans cytoscape + identification des Biotransformations" \
ProteoWizard "Module MSconvert pour convertir les raw dans un format non propriétaire")

if [ $mode = "Metabolites" ]; then
    input_select=$(zenity --file-selection --title="Sélectionner le fichier .clustersummary à traiter")
    ref_select=$(zenity --file-selection --title="Sélectionner le fichier de référence")
    precision_select=$(zenity --forms --title="Sélection de la précision" --text="Le changement de matrice peut avoir un impact sur le RT" --add-entry="Précision de la masse [defaut=0.001]" --add-entry="Précision du RT [defaut=0.1]" --separator=",")
    mass_precision_select=$(echo "$precision_select" | cut -d "," -f1)
    rt_precision_select=$(echo "$precision_select" | cut -d "," -f2)
    ID_metabo -i $input_select -r $ref_select -m $mass_precision_select -t $rt_precision_select
fi

if [ $mode = "Biotransformations" ]; then
    input_select=$(zenity --file-selection --title="Sélectionner le fichier .selfloop à traiter")
    ID_biotransformation -i $input_select
fi

if [ $mode = "InSilico" ]; then
    input_select=$(zenity --file-selection --title="Sélectionner le fichier .txt à traiter")
    option_meta=$(zenity --forms --title="Sélection des options de MetaTrans" --text="add description" --add-entry="Taille min (SMILE) [defaut=5]" --add-entry="Taille max (SMILE) [defaut=120]" --add-entry="Top resultats [defaut=5 : top 10]" --separator=",")
    min=$(echo "$option_meta" | cut -d "," -f1)
    max=$(echo "$option_meta" | cut -d "," -f2)
    beam=$(echo "$option_meta" | cut -d "," -f3)
    type=$(zenity --list --title="Choix du mode de Biotransformers3" --text="Indiquez quel mode utiliser pour Biotransformers3" --column="Mode" --column="Description" \
allHuman Defaut ecbased "add description" cyp450 "add description" phaseII "add description" hgut "add description" superbio "add description" envimicro "add description" )
    predict_metabo -i $input_select -m $min -M $max -t $type -b $beam
fi

if [ $mode = "MetGem" ]; then
    input_select=$(zenity --file-selection --title="Sélectionner le fichier .csv à traiter")
    ID_metgem -i $input_select
fi

if [ $mode = "ProteoWizard" ]; then
    input_select=$(zenity --file-selection --directory --title="Sélectionner le dossier contenant les raw")
    format_select=$(zenity --list --title="Choix du format" --text="Choissisez le format souhaité" --column="Format" --column="Description" \
    mzXML "Defaut" mzML "add description")
    msconvert -i $input_select -f $format_select
fi