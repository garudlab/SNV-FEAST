#!/bin/bash
# 
COUNTER=0

while read -r fam seed uniqueness threshold familythreshold strain minreads study;
do
	if [[ "$study" == *"Suez"* ]]; then
		echo $strain;
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/Private*csv




	    echo /u/scratch/b/briscoel/SuezMontassierFiles/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"*
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
	    cat /u/scratch/b/briscoel/SuezMontassierFiles/data/snps/"$strain"/PrivateSNPs/Private*fam_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*minreads_$minreads*binned_freq.csv > /u/scratch/b/briscoel/SuezMontassierFiles/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv
	    #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/*fam*_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*_freq.csv > /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
	    cat /u/scratch/b/briscoel/SuezMontassierFiles/data/snps/"$strain"/PrivateSNPs/Private*fam_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*minreads_$minreads*counts.csv > /u/scratch/b/briscoel/SuezMontassierFiles/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv

	
	else
		echo $strain;
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/Private*csv


	    echo /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"*
	    #rm -rf /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
	    cat /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/Private*fam_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*minreads_$minreads*binned_freq.csv > /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_binned_freq.csv
	    #cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/*fam*_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*_freq.csv > /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_freq.csv
	    cat /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/Private*fam_$fam*seed_$seed*uniqueness_$uniqueness*threshold_$threshold*familythresh_$familythreshold*minreads_$minreads*counts.csv > /u/home/b/briscoel/project-ngarud/FEASTX/"$study"Files/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_minreads_"$minreads"_counts.csv

	fi

    
    #bzip2  /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv
done < inputs/data_$SGE_TASK_ID.in


#cat /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/PrivateSNPs/*fam*1*seed*"$seed"*uniqueness*"$uniqueness"*threshold*"$threshold"*familythresh*"$familythreshold"*binned_freq.csv > /u/project/ngarud/ngarud/mother-infant/data/snps/"$strain"/CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_binned_freq.csv

#CATTED_fam_"$fam"_seed_"$seed"_uniqueness_"$uniqueness"_threshold_"$threshold"_familythresh_"$familythreshold"_counts.csv