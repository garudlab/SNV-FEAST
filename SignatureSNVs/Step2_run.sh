#!/bin/bash
. /u/local/Modules/default/init/modules.sh
#source /u/home/b/briscoel/project-ngarud/miniconda2/bin/activate /u/home/b/briscoel/project-ngarud/miniconda2/envs/python3
#source /u/home/b/briscoel/project-ngarud/miniconda2/bin/activate /u/home/b/briscoel/project-ngarud/miniconda2/envs/python3
module load anaconda3/2020.11
while read -r arg_1 arg_2 arg_3 arg_4 arg_5 arg_6 arg_7 arg_8 arg_9; do
	python Step2_get_sigs.py --seed $arg_1 --uniqueness $arg_2 --threshold $arg_3 --family_threshold $arg_4 --strain $arg_5 --min_reads $arg_6 --start_index $arg_7 --end_index $arg_8 --study $arg_9
done < inputs/data_$SGE_TASK_ID.in






