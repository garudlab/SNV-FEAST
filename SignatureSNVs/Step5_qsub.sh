
# ./Step4_qsub.sh ../high_prev_strains.txt ../../BackhedFiles/Pipe_CompleteFams.csv 5000 2 1 0.8 0.1
#./Step5_qsub.sh ../../BackhedFiles/Pipe_CompleteFams.csv 1 3 1 0.8 0.1
# ./Step5_qsub.sh ../../BackhedFiles/Pipe_CompleteFams.csv 100 4 1 0.5 0.1

first_count_input=$2
seed=$3
uniqueness=$4
threshold=$5
familythreshold=$6
minreads=$7
study=$8
COUNTER="$(($first_count_input-1))"
while IFS= read -r fam;
do
	COUNTER=$((COUNTER + 1));        
	echo $strain;
	echo $fam
	echo "$fam $seed $uniqueness $threshold $familythreshold $minreads $study" > inputs/data_$COUNTER.in; 
done < $1


if [ "$#" -ge 9 ]; then
    echo "9 parames"
    qsub -hold_jid $9 -cwd -V -N s"$seed"step5 -e misc -o misc -l h_data=6G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step5_zip.sh"

    # 12G




else 
	echo "less than 9 param"
	
	qsub -cwd -V -N s"$seed"step5 -e misc -o misc -l h_data=6G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step5_zip.sh"

fi



#qsub -hold_jid 6576791 -cwd -V -N ZipOnly -l h_data=12G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./Step5_zip.sh"


# qsub -cwd -V -N ZipOnly -l h_data=12G,time=24:00:00,highp -b y -t 55:55 "./Step5_zip.sh"

