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

echo "
===============================================================================================================
||                                                                                           
||   Process of ${molecule} : ${smile} by Biotransformers3
||                                                                                            
===============================================================================================================
"

#Singularity version
singularity exec https://depot.galaxyproject.org/singularity/biotransformer%3A3.0.20230403--hdfd78af_0 biotransformer \
-b "${type}" \
-k "pred" \
-cm 3 \
-s "${nstep}" \
-ismi "${smile}" \
-ocsv "${DirOutput}${molecule}_Biotransformer3.csv"

#-a : annotate Search PuChem for each product, and annotate with CID and synonyms

#Download version
#java -jar BioTransformer3.0_20230525.jar -b ${type} -k pred -cm 3 -ismi "${smile}" -ocsv "${DirOutput}${molecule}_Biotransformer3.csv"

#Changement csv format
csvformat -D ";" "${DirOutput}${molecule}_Biotransformer3.csv" | gawk -v RS='"' 'NR % 2 == 0 { gsub(/\n/, "") } { printf("%s%s", $0, RT) }' > "${DirOutput}${molecule}_Biotransformer3_brut.csv"

#Remove useless columns
molecule_nb=0
while read line
do
    if [[ $line =~ $pat_bioTrans ]]; then
        :
    else
        ((molecule_nb+=1))
        info=$(echo $line | cut -d\; -f3,6,8,14,16,17)
        echo "Molecule_${molecule_nb};${info}" >> ${DirOutput}${molecule}_Biotransformer3_brut2.csv
    fi
done < ${DirOutput}${molecule}_Biotransformer3_brut.csv
sed 's/;/,/g' ${DirOutput}${molecule}_Biotransformer3_brut2.csv > ${DirOutput}${molecule}_Biotransformer3_brut3.csv
rm ${DirOutput}${molecule}_Biotransformer3_brut.csv ${DirOutput}${molecule}_Biotransformer3_brut2.csv ${DirOutput}${molecule}_Biotransformer3.csv

###Structure construction
python3 $Script_SmitoStr -i "${DirOutput}${molecule}_Biotransformer3_brut3.csv"
mkdir ${DirFigBioTrans}/${molecule}
mv Molecule*.jpeg ${DirFigBioTrans}/${molecule}

#Calcul Mass+H from brut formula
while read line
do
    formula=$(echo $line | cut -d\, -f3)
    python3 $Script_massFromFormula $formula
done < "${DirOutput}${molecule}_Biotransformer3_brut3.csv" >> "${DirOutput}${molecule}_Biotransformer3_mass.csv"

#Merging columns
paste -d ',' ${DirOutput}${molecule}_Biotransformer3_brut3.csv ${DirOutput}${molecule}_Biotransformer3_mass.csv > ${DirOutput}${molecule}_Biotransformer3_VF.csv
rm ${DirOutput}${molecule}_Biotransformer3_brut3.csv ${DirOutput}${molecule}_Biotransformer3_mass.csv

#Creating entetes
echo "Métabolites,SMILES,FormuleBrute,ALogP,Pathway,Enzyme,Systeme,Masse(+H)" > ${DirOutput}${molecule}_BioTransformer3.csv
while read line
do 
    echo $line
done < "${DirOutput}${molecule}_Biotransformer3_VF.csv" >> ${DirOutput}${molecule}_BioTransformer3.csv
rm ${DirOutput}${molecule}_Biotransformer3_VF.csv

conda deactivate
}