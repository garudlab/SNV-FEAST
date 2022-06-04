
# ./Step4_qsub.sh ../high_prev_strains.txt ../../BackhedFiles/Pipe_CompleteFams.csv 5000 2 1 0.8 0.1
#./Step4_qsub.sh ../high_prev_strains.txt ../../BackhedFiles/Pipe_CompleteFams.csv 3000 3 1 0.8 0.1 

first_count_input=$3
seed=$4
uniqueness=$5
threshold=$6
familythreshold=$7
COUNTER="$(($first_count_input-1))"
minreads=$8
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
    qsub -hold_jid ${10} -cwd -V -N Seed"$seed"step4 -e mis -o mis -l h_data=8G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./Step4_catzip_acrossstrain.sh"



else 
	echo "less than 10 param"
	
	qsub -cwd -V -N Seed"$seed"s4 -e mis -o mis -l h_data=8G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./Step4_catzip_acrossstrain.sh"

	#qsub -cwd -V -N Seed6s4 -e misc -o misc -l h_data=8G,time=24:00:00,highp -b y -t 7365-7999 "./Step4_catzip_acrossstrain.sh"

fi

#qsub -hold_jid 6576762 -cwd -V -N CatAcrosszip -l h_data=8G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step4_catzip_acrossstrain.sh"

# qsub -cwd -V -N CatAcrosszip -l h_data=8G,time=24:00:00 -b y -t 1033:1843 "./Step4_catzip_acrossstrain.sh"

# qsub -cwd -V -N CatAcrosszip -l h_data=8G,time=24:00:00 -b y -t 1:1 "./Step4_catzip_acrossstrain.sh"