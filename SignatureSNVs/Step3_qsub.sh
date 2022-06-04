#!/bin/bash
# example ./qsub_cat_and_zip.sh strain_files.txt 100 0.5
# ./Step3_qsub.sh ../test_strain_files.txt ../../BackhedFiles/Pipe_CompleteFams.csv 1 2 1 0.8 0.1
# ./Step3_qsub.sh ../high_prev_strains.txt ../../BackhedFiles/Pipe_CompleteFams.csv 1200 2 1 0.8 0.1
# ./Step3_qsub.sh ../high_prev_strains.txt ../../BackhedFiles/Pipe_CompleteFams.csv 1 0 1 0.8 0.1 10

first_count_input=$3
seed=$4
uniqueness=$5
threshold=$6
familythreshold=$7
minreads=$8
COUNTER="$(($first_count_input-1))"
study=$9



while IFS= read -r strain;
do
	while IFS= read -r fam;
	do
		COUNTER=$((COUNTER + 1));        
		echo $strain;
		echo $fam
		echo "$fam $seed $uniqueness $threshold $familythreshold $strain $minreads $study" > inputs/data_$COUNTER.in; 
	done < $2
done < $1

if [ "$#" -ge 10 ]; then
    echo "9 parames"
    qsub  -hold_jid ${10} -cwd -V -N S"$seed"step3 -e mis -o mis -l h_data=6G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./Step3_cat_withinstrain.sh"


else 
	echo "less than 9 param"
	qsub -cwd -V -N S"$seed"s3 -e mis -o mis -l h_data=6G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step3_cat_withinstrain.sh"
fi
#qsub -hold_jid 6576728 -cwd -V -N Catzip -l h_data=6G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step3_cat_withinstrain.sh"
