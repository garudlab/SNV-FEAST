#!/bin/bash

# seed=$3
# uniqueness=$4
# threshold=$5
# familythreshold=$6

while read -r fam seed uniqueness threshold familythreshold strain minreads study;
do
    if [[ "$study" == *"Suez"* ]]; then
        echo $strain;

        echo CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_binned_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_counts.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_counts.csv
        
        echo "not 0 reads"
        cat /u/scratch/b/briscoel/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv >> /u/scratch/b/briscoel/"$study"Files/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
        cat /u/scratch/b/briscoel/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv >> /u/scratch/b/briscoel/"$study"Files/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv
    else
        echo $strain;

        echo CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_binned_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_counts.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_counts.csv
        
        echo "not 0 reads"
        cat /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv
        #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
        cat /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv
 
    fi


    #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_$fam_seed_$seed_uniqueness_$uniqueness_threshold_$threshold_familythresh_$familythreshold_binned_freq.csv >> /u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
    #/u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_counts.csv
    #bzip2  /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv

done < inputs/data_$SGE_TASK_ID.in


#qsub -cwd -V -N CatzipBIG -l h_data=16G,time=24:00:00,highp -b y "./Step4_catzip_acrossstrain.sh"
