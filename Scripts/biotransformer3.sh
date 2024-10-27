#!  /usr/bin/bash

biotransformer_function () {

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
        -t|--type)
            type="$2"
            shift
            shift
            ;;
        -b|--biodir)
            DirBiotrans="$2"
            shift
            shift
            ;;
        -r|--script)
            DirScripts="$2"
            shift
            shift
            ;;
        -n|--nstep)
            nstep="$2"
            shift
            shift
            ;;
        -c|--cmode)
            cmode="$2"
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

DirFigBioTrans=${DirOutput}BioTrans_Figures
mkdir $DirFigBioTrans

#cd ${DirBiotrans}

pat_bioTrans="^InChI;InChIKey"

for indice in ${!tab_molecule[@]}
do
    echo "
 ===============================================================================================================
||                                                                                           
||   Process of ${tab_molecule[${indice}]} : ${tab_smiles[${indice}]} by Biotransformers3
||                                                                                            
 ===============================================================================================================
    "
    
    #Singularity version
    singularity exec https://depot.galaxyproject.org/singularity/biotransformer%3A3.0.20230403--hdfd78af_0 biotransformer \
    -b "${type}" \
    -k "pred" \
    -cm 3 \
    -s "${nstep}" \
    -ismi "${tab_smiles[${indice}]}" \
    -ocsv "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3.csv"

    #-a : annotate Search PuChem for each product, and annotate with CID and synonyms

    #Download version
    #java -jar BioTransformer3.0_20230525.jar -b ${type} -k pred -cm 3 -ismi "${tab_smiles[${indice}]}" -ocsv "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3.csv"

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

    ###Structure construction
    python3 $Script_SmitoStr -i "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv"
    mkdir ${DirFigBioTrans}/${tab_molecule[${indice}]}
    mv Molecule*.jpeg ${DirFigBioTrans}/${tab_molecule[${indice}]}

    #Calcul Mass+H from brut formula
    while read line
    do
        formula=$(echo $line | cut -d\, -f3)
        python3 $Script_massFromFormula $formula
    done < "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv" >> "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv"

    #Merging columns
    paste -d ',' ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv > ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_brut3.csv ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_mass.csv

    #Creating entetes
    echo "Métabolites,SMILES,FormuleBrute,ALogP,Pathway,Enzyme,Systeme,Masse(+H)" > ${DirOutput}${tab_molecule[${indice}]}_BioTransformer3.csv
    while read line
    do 
        echo $line
    done < "${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv" >> ${DirOutput}${tab_molecule[${indice}]}_BioTransformer3.csv
    rm ${DirOutput}${tab_molecule[${indice}]}_Biotransformer3_VF.csv

done

conda deactivate
}
