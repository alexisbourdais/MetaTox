#!  /usr/bin/bash

sygma_function () {

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
        -m|--molecule)
            molecule="$2"
            shift
            shift
            ;;
        -s|--smile)
            smile="$2"
            shift
            shift
            ;;
        -o|--outdir)
            DirOutput="$2"
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

conda activate my-rdkit-env

DirFigSygma=${DirOutput}Sygma_Figures
mkdir $DirFigSygma

pat_pathway="^>  <Pathway>"
pat_score="^>  <Score>"

echo "
#####   Process of $molecule : $smile by Sygma   #####
"

singularity run docker://3dechem/sygma ${smile} >> "${DirOutput}${molecule}_Sygma.sdf"

#--phase1 Number of phase 1 cycles (default: 1)
#--phase2 Number of phase 2 cycles (default: 1)
#-o, --outputtype Molecule output type (default: sdf)

###Conversion des sdf en smiles
python3 ${Script_SdftoSmi} ${DirOutput}${molecule}_Sygma.sdf
mv ${PWD}/smiles.txt ${DirOutput}Prediction_sygma_${molecule}_smiles.txt

###Réalisation des structures
mkdir ${DirFigSygma}/${molecule}
cd ${DirFigSygma}/${molecule}
python3 $Script_SmitoStr -i ${DirOutput}Prediction_sygma_${molecule}_smiles.txt

#Ajout des scores & pathway du fichier sdf dans le fichier txt
while read line
do
    if [[ $line_pre == "pathway" ]]; then
        echo $line | sed 's/,/;/g' >> "${DirOutput}Prediction_sygma_${molecule}_Path.txt"
        line_pre=""
    fi

    if [[ $line_pre == "score" ]]; then
        echo $line >> "${DirOutput}Prediction_sygma_${molecule}_Score.txt"
        line_pre=""
    fi

    if [[ $line =~ $pat_pathway ]]; then
        line_pre=pathway
    fi

    if [[ $line =~ $pat_score ]]; then
        line_pre=score
    fi
done < "${DirOutput}${molecule}_Sygma.sdf"

paste -d ',' ${DirOutput}Prediction_sygma_${molecule}_smiles.txt ${DirOutput}Prediction_sygma_${molecule}_Path.txt ${DirOutput}Prediction_sygma_${molecule}_Score.txt > ${DirOutput}${molecule}_Sygma_SmilesPathScore.txt
rm ${DirOutput}Prediction_sygma_${molecule}_smiles.txt ${DirOutput}Prediction_sygma_${molecule}_Path.txt ${DirOutput}Prediction_sygma_${molecule}_Score.txt

#Convert smiles into brut formula
while read line
do
    smiles=$(echo $line | cut -d\, -f2)
    python3 $Script_FormulaFromSmiles $smiles
done < ${DirOutput}${molecule}_Sygma_SmilesPathScore.txt >> ${DirOutput}${molecule}_Sygma_formuleBrute.txt

#Calcul Mass+H from brut formula
while read line
do
    formula=$(echo $line)
    python3 $Script_massFromFormula $formula
done < ${DirOutput}${molecule}_Sygma_formuleBrute.txt >> ${DirOutput}${molecule}_Sygma_mass.txt

#Compilation des colonnes
paste -d ',' ${DirOutput}${molecule}_Sygma_SmilesPathScore.txt ${DirOutput}${molecule}_Sygma_formuleBrute.txt ${DirOutput}${molecule}_Sygma_mass.txt > ${DirOutput}${molecule}_SygmaVF.txt
rm ${DirOutput}${molecule}_Sygma_SmilesPathScore.txt ${DirOutput}${molecule}_Sygma_formuleBrute.txt ${DirOutput}${molecule}_Sygma_mass.txt 

#Création des entêtes
echo "Métabolites,SMILES,Pathway,Score,FormuleBrute,Masse(+H)" > ${DirOutput}${molecule}_Sygma.csv
while read line
do 
    echo $line | sed "s/;//g"
done < ${DirOutput}${molecule}_SygmaVF.txt >> ${DirOutput}${molecule}_Sygma.csv
rm ${DirOutput}${molecule}_SygmaVF.txt
}