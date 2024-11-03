#!  /usr/bin/bash

#########
# TO DO #
#########

#boucle while molecule molecule pour chaque ${nom molecule}_metapred-sygma-metatrans.csv

############################
### Barre de progression ###
############################

bar_size=40
bar_char_done="#"
bar_char_todo="-"
bar_percentage_scale=2

function show_progress {
    current="$1"
    total="$2"

    # calculate the progress in percentage 
    percent=$(bc <<< "scale=$bar_percentage_scale; 100 * $current / $total" )
    # The number of done and todo characters
    done=$(bc <<< "scale=0; $bar_size * $percent / 100" )
    todo=$(bc <<< "scale=0; $bar_size - $done" )

    # build the done and todo sub-bars
    done_sub_bar=$(printf "%${done}s" | tr " " "${bar_char_done}")
    todo_sub_bar=$(printf "%${todo}s" | tr " " "${bar_char_todo}")

    # output the bar
    echo -ne "\rProgress : [${done_sub_bar}${todo_sub_bar}] ${percent}%"

    if [ $total -eq $current ]; then
        echo -e "\nDONE"
    fi
}

############################
###     Fonction index   ###
############################

function get_index() {
for i in "${!smile_tab[@]}"; do
   if [[ "${smile_tab[$i]}" = "${1}" ]]; then
       echo "${i}";
   fi
done
}

#####################
###     Program   ###
#####################

dir_res="Results_Prediction"
molecule="Nicotine"

sygma_file="$dir_res/${molecule}_Sygma.csv"
trans_file="$dir_res/${molecule}_Metatrans.csv"
pred_file="$dir_res/${molecule}_metapred.csv"

pattern="^Métabolites*"
results_file="Fusion_results.csv"

# Sygma table
smile_tab=()
sygma_pathway_tab=()
sygma_score_tab=()
formuleBrute_tab=()
mass_tab=()
metapred_tab=()
metatrans_tab=()
sygma_tab=()

#############
### Sygma ###
#############

tasks_in_total=$( nl "$sygma_file" | wc -l )
current_task=0

while read -r line
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total

    if [[ $line =~ $pattern ]]; then
        :

    else
        id=$(echo "$line" | cut -d',' -f1)
        smile=$(echo "$line" | cut -d',' -f2)
        pathway=$(echo "$line" | cut -d',' -f3)
        score=$(echo "$line" | cut -d',' -f4)
        formulebrute=$(echo "$line" | cut -d',' -f5)
        masse=$(echo "$line" | cut -d',' -f6)

        smile_tab[${#smile_tab[@]}]=${smile}
        sygma_pathway_tab[${#sygma_pathway_tab[@]}]=${pathway}
        sygma_score_tab[${#sygma_score_tab[@]}]=${score}
        formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
        mass_tab[${#mass_tab[@]}]=${masse}
        index=$(get_index $smile)
        sygma_tab[${index}]="+"

    fi
done < "$sygma_file"

#################
### MetaTrans ###
#################

tasks_in_total=$( nl "$trans_file" | wc -l )
current_task=0

while read -r line
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total

    if [[ $line =~ $pattern ]]; then
        :

    else
        id=$(echo "$line" | cut -d',' -f1)
        smile=$(echo "$line" | cut -d',' -f2)
        formulebrute=$(echo "$line" | cut -d',' -f3)
        masse=$(echo "$line" | cut -d',' -f4)

        #si la sequence deja presente dans le tableau
        if [[ $(echo ${smile_tab[@]} | fgrep -w $smile) ]]; then
            index=$(get_index $smile)
            metatrans_tab[${index}]="+"

        #si la sequence n'est pas encore dans le tableau, l'ajoute
        else
            smile_tab[${#smile_tab[@]}]=${smile}
            formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
            mass_tab[${#mass_tab[@]}]=${masse}
            index=$(get_index $smile)
            metatrans_tab[${index}]="+"
        fi
    fi

done < "$trans_file"

######################
### Meta-Predictor ###
######################

tasks_in_total=$( nl "$pred_file" | wc -l )
current_task=0

while read -r line
do
    ((current_task+=1))
    show_progress $current_task $tasks_in_total

    if [[ $line =~ $pattern ]]; then
        :

    else
        id=$(echo "$line" | cut -d',' -f1)
        smile=$(echo "$line" | cut -d',' -f2)
        formulebrute=$(echo "$line" | cut -d',' -f3)
        masse=$(echo "$line" | cut -d',' -f4)

        #si la sequence deja presente dans le tableau
        if [[ $(echo ${smile_tab[@]} | fgrep -w $smile) ]]; then
            index=$(get_index $smile)
            metapred_tab[${index}]="+"

        #si la sequence n'est pas encore dans le tableau, l'ajoute
        else
            smile_tab[${#smile_tab[@]}]=${smile}
            formuleBrute_tab[${#formuleBrute_tab[@]}]=${formulebrute}
            mass_tab[${#mass_tab[@]}]=${masse}
            index=$(get_index $smile)
            metapred_tab[${index}]="+"
        fi
    fi
done < "$pred_file"

###################
### Comparaison ###
###################

# Entete
echo -e "FormuleBrute\tMasse(+H)\tSmile\tSygma_pathway\tSygma_score\tSygma\tMetaTrans\tMetaPredictor" > ${results_file}

for indice in ${!smile_tab[@]}
do
    smile="${smile_tab[${indice}]}"
    formulebrute="${formuleBrute_tab[${indice}]}"
    masse="${mass_tab[${indice}]}"
    score="${sygma_score_tab[${indice}]}"
    pathway="${sygma_pathway_tab[${indice}]}"
    metapred="${metapred_tab[${indice}]}"
    metatrans="${metatrans_tab[${indice}]}"
    sygma="${sygma_tab[${indice}]}"

    echo -e "${formulebrute}\t${masse}\t${smile}\t${pathway}\t${score}\t${sygma}\t${metatrans}\t${metapred}" >> ${results_file}

done