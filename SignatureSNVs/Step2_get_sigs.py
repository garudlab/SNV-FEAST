# FUNCTION: GET SIGNATURE SNPS family by family

# EXAMPLE
# python Step1.py --local_bool 0 --tester_mode 1 --uniqueness 1 --threshold 0.20 --family_threshold 0.5 --strain Bacteroides --min_reads 20 --start_index 1 --end_index 5
# python Step1.py --local_bool --tester_mode --uniqueness 1 --threshold 0.20 --family_threshold 0.5 --strain Bacteroides --min_reads 20 --start_index 1 --end_index 5
# python Step2_get_sigs.py --seed 0 --uniqueness 1 --threshold 0.1 --family_threshold 0.1 --strain Escherichia_coli_53492 --min_reads 20 --start_index 30000 --end_index 60000

# python Step2_get_sigs.py --seed 4 --uniqueness 1 --threshold 0.8 --family_threshold 0.1 --strain Bacteroides_uniformis_57318 --min_reads 10 --start_index 1 --end_index 200000
# started at 4 PM Jan 14
# END EXAMPLE

#python Step2_get_sigs.py --seed 1 --uniqueness 1 --threshold 1.0 --family_threshold 0.1 --strain Klebsiella_pneumoniae_54788 --min_reads 10 --start_index 1 --end_index 20000 --study Brooks   


#python Step2_get_sigs.py --seed 3 --uniqueness 1 --threshold 1 --family_threshold 0.1 --strain Bacteroides_vulgatus_57955 --min_reads 5 --start_index 1 --end_index 2000 --study SimulationB   

# python Step2_get_sigs.py --seed 2 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Bacteroides_vulgatus_57955 --min_reads 5 --start_index 1 --end_index 200000 --study SimulationE

# python Step2_get_sigs.py --seed 1 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Pelagibacter_ubique_58575 --min_reads 5 --start_index 1 --end_index 2000 --study Tara  


# python Step2_get_sigs.py --seed 17 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Bacteroides_vulgatus_57955  --min_reads 10 --start_index 200001 --end_index 220000 --study Backhed

#python Step2_get_sigs.py --seed 2 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Bacteroides_uniformis_57318  --min_reads 10 --start_index 1 --end_index 30 --study Dummy  

#python Step2_get_sigs.py --seed 21 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Roseburia_inulinivorans_61943 --min_reads 10 --start_index 1 --end_index 2000 --study SimulationD  
# python Step2_get_sigs.py --seed 4 --uniqueness 1 --threshold 1 --family_threshold 1.0 --strain Bacteroides_vulgatus_57955 --min_reads 10 --start_index 1 --end_index 200 --study SuezMontassier


# needed moduels
import argparse, sys
import numpy as np
import pandas as pd
import os
from helpers import make_np_array_from_file, get_private_snps, score_snps
import time
from itertools import compress
import random
#paramsbzcat f
print("NEW CODE")
parser=argparse.ArgumentParser()

parser.add_argument('--local_bool', help='Running locally (0 or 1)',default=False, action='store_true')
parser.add_argument('--tester_mode', help='Minimal tester mode? (0 or 1)',default=False, action='store_true')
parser.add_argument('--seed', help='Seed for random fam draw',type=int)
parser.add_argument('--uniqueness', help='How many mothers can have that SNP?',type=int)
parser.add_argument('--threshold', help='What is the minimum allele frequency',type =float)
parser.add_argument('--family_threshold', help='What is the minimum allele frequency for non mothers', type = float)
parser.add_argument('--strain', help='Strain for this job')
parser.add_argument('--min_reads', help='Minimum read depth for a valid allele frequency',type = int)
parser.add_argument('--start_index', help='Starting index for checking the strain file because big file', type = int)
parser.add_argument('--end_index', help='End index', type = int)
parser.add_argument('--study', help='Study to study',type=str)


args=parser.parse_args()
print(args)

get_summary_stats = True

run_clean = False # delete files after
local_bool = args.local_bool
tester_mode = args.tester_mode
seed = args.seed
uniqueness = args.uniqueness
threshold = args.threshold
familythreshold = args.family_threshold
strain = args.strain
min_reads = args.min_reads # minimum read depth
start_index = args.start_index #following python system 0 based
end_index = args.end_index

study = args.study


print("uniqueness" + str(uniqueness))
print("thresh" + str(threshold))
custom_snv = False
if local_bool:
	input_dir =  '/Users/leahbriscoe/Documents/FEASTX/BackhedFiles/'
	snp_dir = input_dir + 'snps_' + strain + '/'

else:

	if study == "Backhed":
		input_dir =  '/u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/'
		snp_dir = '/u/project/ngarud/daisyche/mother_infant/data/snps/' + strain + '/'
		out_snp_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/BackhedFiles/data/snps/' + strain 
	elif study == "Brooks":
		input_dir =  '/u/home/b/briscoel/project-ngarud/FEASTX/BrooksFiles/'
		snp_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/' + study + 'Files/midas_output_snps/' + strain + '/'
		#'/u/home/b/briscoel/project-ngarud/NICU_data/midas_output_snps/' + strain + '/'
		snp2_dir = '/u/project/ngarud/daisyche/mother_infant/data/snps/' + strain + '/'

		out_snp_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/BrooksFiles/data/snps/' + strain 
	elif study == "SuezMontassier":
		if seed == 4:
			custom_snv = True
		fst_dir = "/u/scratch/b/briscoel/SuezMontassierFiles/Fst_data/"
		input_dir =  '/u/scratch/b/briscoel/SuezMontassierFiles/'
		snp_dir = '/u/scratch/b/briscoel/' + study + 'Files/midas_output_snps/' + strain + '/'
		#'/u/home/b/briscoel/project-ngarud/NICU_data/midas_output_snps/' + strain + '/'
		snp2_dir = '/u/scratch/b/briscoel/' + study + 'Files/midas_output_snps/'  + strain + '/'

		out_snp_dir = '/u/scratch/b/briscoel/SuezMontassierFiles/data/snps/' + strain 
	else:
		input_dir =  '/u/home/b/briscoel/project-ngarud/FEASTX/' + study + 'Files/'
		snp_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/' + study + 'Files/midas_output_snps/' + strain + '/'
		out_snp_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/' + study + 'Files/data/snps/' + strain 
		

	if not os.path.exists(out_snp_dir ):
		os.makedirs(out_snp_dir )
	out_snp_dir = out_snp_dir + '/'
	
	

	#snp_dir = '/u/project/ngarud/ngarud/mother-infant/data/snps/' + strain + '/'


## METADATA
family_record = pd.read_csv(input_dir +"Pipe_ChosenFamily_seed" + str(seed)  + ".csv")
# if study == "SuezMontassier":
# 	family_record = pd.read_csv(input_dir +"Pipe_ChosenFamily_seed" + str(seed)  + ".csv",skip_blank_lines=True)

print("Family record")
print(family_record)
print(family_record.loc[0])

all_families_samples = family_record.iloc[:,1:(family_record.shape[1])].values
all_families_samples = all_families_samples.flatten()

# if study == "Backhed":
# 	all_mom_samples = family_record.iloc[:,4:(family_record.shape[1])].values
if study =="Li":
	if seed%2 == 1:
		all_mom_samples = family_record.iloc[:,6:(family_record.shape[1])].values
	else:
		all_mom_samples = family_record.iloc[:,5:(family_record.shape[1])].values
else:
	all_mom_samples = family_record.iloc[:,2:(family_record.shape[1])].values
	all_infant_samples = family_record.iloc[:,1:2].values

all_mom_samples = sorted(list(set(all_mom_samples.flatten())))

all_infant_samples = sorted(list(set(all_infant_samples.flatten())))

print(all_families_samples)
print("all mom samples")
print(all_mom_samples)
print("All infant samples")
print(all_infant_samples )

family_record.index = family_record['family_id']
## DATA PREP

file_string = "PrivateSNPs_ind_" + str(start_index) + "_" + str(end_index) + "_uniqueness_" + str(uniqueness) +\
"_threshold_" + str(threshold) + "_familythresh_" + str(familythreshold)

all_file_string = "PrivateSNPs" + "_seed_" + str(seed) + "_ind_" + str(start_index) + "_" + str(end_index) + "_uniqueness_" + str(uniqueness) +\
"_threshold_" + str(threshold) + "_familythresh_" + str(familythreshold) + "_minreads_" + str(min_reads)
depth_path = snp_dir + "snps_depth.txt.bz2"

freq_path = snp_dir + "snps_ref_freq.txt.bz2"



if study == "Brooks": # because comning NICU samples with baby samples from another folder
	print("Brooks")
	depth_path2 = snp2_dir + "snps_depth.txt.bz2"
	freq_path2 = snp2_dir + "snps_ref_freq.txt.bz2"


	snp_depth_matrix2,depth_row_snp_names2,depth_col_sample_names2 = make_np_array_from_file(depth_path2,start_index,end_index,all_samples=all_families_samples)
	snp_freq_matrix2,freq_row_snp_names2,freq_col_sample_names2 = make_np_array_from_file(freq_path2,start_index,end_index,all_samples=all_families_samples)
	snp_freq_matrix2 = 1- snp_freq_matrix2 # to make it alternative allele frequency

	#if snp_freq
	snp_depth_matrix1,depth_row_snp_names1,depth_col_sample_names1 = make_np_array_from_file(depth_path,start_index,end_index,all_samples=all_families_samples)
	snp_freq_matrix1,freq_row_snp_names1,freq_col_sample_names1 = make_np_array_from_file(freq_path,start_index,end_index,all_samples=all_families_samples)
	snp_freq_matrix1 = 1- snp_freq_matrix1 # to make it alternative allele frequency
	
	snp_freq_matrix = np.concatenate((snp_freq_matrix2, snp_freq_matrix1), axis=1)
	snp_depth_matrix = np.concatenate((snp_depth_matrix2, snp_depth_matrix1), axis=1)
	
	depth_col_sample_names = depth_col_sample_names2 + depth_col_sample_names1
	freq_row_snp_names = [strain + "|" + i for i in freq_row_snp_names2]
	freq_col_sample_names = freq_col_sample_names2 + freq_col_sample_names1

	depth_row_snp_names = depth_row_snp_names1
	snp_depth_matrix2 = []
	snp_freq_matrix2 = []
	snp_depth_matrix1 = []
	snp_freq_matrix1  = []



else:
	if custom_snv:
		fst_stats = pd.read_csv(fst_dir + "VEC-Fst-" + strain +  "-Regular.txt",index_col=0,sep="\t",squeeze=True)
		print("fst_stats")
		print(fst_stats)
		print("custom features")
		print(fst_stats[fst_stats>=0.25].index)
		custom_feat = fst_stats[fst_stats>=0.25].index
		snp_depth_matrix,depth_row_snp_names,depth_col_sample_names = make_np_array_from_file(depth_path,start_index,end_index,all_samples=all_families_samples,custom_features = custom_feat)
		snp_freq_matrix,freq_row_snp_names,freq_col_sample_names = make_np_array_from_file(freq_path,start_index,end_index,all_samples=all_families_samples,custom_features = custom_feat)
	
	else:
		snp_depth_matrix,depth_row_snp_names,depth_col_sample_names = make_np_array_from_file(depth_path,start_index,end_index,all_samples=all_families_samples)
		snp_freq_matrix,freq_row_snp_names,freq_col_sample_names = make_np_array_from_file(freq_path,start_index,end_index,all_samples=all_families_samples)
	snp_freq_matrix = 1- snp_freq_matrix # to make it alternative allele frequency

	freq_row_snp_names = [strain + "|" + i for i in freq_row_snp_names]

infant_bit_mask = [s in all_infant_samples for s in freq_col_sample_names]
snp_freq_matrix_infants = snp_freq_matrix[:,infant_bit_mask]
print(snp_freq_matrix_infants.shape)
snp_presence = snp_freq_matrix_infants  > 0 
snp_prevalence = snp_presence.sum(axis=1)/snp_presence.shape[1]

print("Print prevalene")
print(snp_prevalence.shape)
print(snp_prevalence[0:100])

# snp_prevalence < threshold means maxing out the prevalence(hopw common snps are among infants at no more than 1)
# basicallt we want to take only snps that are observed at least in one of the infants
# and we still want snips that are not in mother at all, so really the mother prevalence is irrelevant (if we want to approximate uniquenss of infant properly and get unknowns)
# previously np.array(snp_prevalence < threshold) & np.array(snp_prevalence > 0)

prevalence_filter = np.array(snp_prevalence <= threshold) & np.array(snp_prevalence > 0)
print("sum filter")
print(sum(prevalence_filter))
print(prevalence_filter)


freq_row_snp_names = np.array(freq_row_snp_names)[prevalence_filter]
depth_row_snp_names = np.array(depth_row_snp_names)[prevalence_filter]


snp_freq_matrix = snp_freq_matrix[prevalence_filter,:]
snp_depth_matrix = snp_depth_matrix[prevalence_filter,:]


#feature_filter = min_depth_bit_mask.sum(axis=1) == min([2,num_cols] ) # if there are at least 2 samples, at least both should have that minimum count of observation

print("shape after filter")
print(snp_depth_matrix.shape)
print(snp_depth_matrix[1:10])
 
	


## PSUEDOCODE
# for each child
# 	draw 3 random moms
# 	get mothers private snps

if study == "Backhed":
	unique_family_ids = [int(f) for f in family_record['family_id']]
if study == "SimulationT":
	if seed %2 == 0:
		species_presence_key = pd.read_csv(input_dir +"transmission_offhand.csv",index_col = 0)
	else:
		species_presence_key = pd.read_csv(input_dir +"transmission_ncm.csv",index_col = 0)
# if study == "SimulationD":
# 	if seed  >= 8:
# 		species_presence_key = pd.read_csv(input_dir +"transmission_offhand.csv",index_col = 0)
	
unique_family_ids= [f for f in family_record['family_id']]

complete_families = []

all_fam_LR = dict()
all_fam_present = dict()
fam_count = 0

print("FAM record")
print(family_record)
print("UNIQUE_FAMILH ID")
print(unique_family_ids)

for fam in unique_family_ids:
	fam_count += 1
	print("FAM")
	print(fam)


	## CHECK if species is evaluated for this infant?
	if study == "SimulationT":  #or ((study == "SimulationD") & (seed  >= 8)):
		print("Transmission record")
		print(species_presence_key.loc[fam][strain])
		if species_presence_key.loc[fam][strain] < 1:
			print("skip " + str(fam))
			continue
	
	start_time = time.time()
	file_string = "PrivateSNPs_fam_" + str(fam) + "_seed_" + str(seed) + "_ind_" + str(start_index) + "_" + str(end_index) + "_uniqueness_" + str(uniqueness) +\
"_threshold_" + str(threshold) + "_familythresh_" + str(familythreshold) + "_minreads_" + str(min_reads)


	all_fam = list(family_record.loc[fam][1:family_record.shape[1]])
	#all_fam = list(family_record.loc[fam].iloc[:,1:family_record.shape[1]])
	# print("all_fam")
	# print(all_fam)
	# print("all_fam length")
	# print(len(all_fam))

	# if study == "Backhed":
	# 	all_moms = list(family_record.loc[fam][4:family_record.shape[1]])
	if study == "Li":
		if seed %2 == 1:
			all_moms = list(family_record.loc[fam][6:family_record.shape[1]])
		else:
			all_moms = list(family_record.loc[fam][5:family_record.shape[1]])
	else:
		all_moms = list(family_record.loc[fam][2:family_record.shape[1]])



	print("all moms")
	print(all_moms)
	print("allfam")
	print(all_fam)
	print("freq_col_sample_names")
	print(freq_col_sample_names)

	mother_bit_mask = [s in all_moms for s in freq_col_sample_names]
	family_bit_mask = [s in all_fam for s in freq_col_sample_names]
	familywise_bit_mask = [s in freq_col_sample_names for s in all_fam]
	all_fam_present[fam]  = pd.Series(familywise_bit_mask )
	# if study == "Backhed":
	# 	infant_bit_mask = [s in all_fam[0:3] for s in freq_col_sample_names]
	if study == "Li":
		if seed %2 == 1:
			infant_bit_mask = [s in all_fam[0:5] for s in freq_col_sample_names]
		else:
			infant_bit_mask = [s in all_fam[0:4] for s in freq_col_sample_names]

	else:
		infant_bit_mask = [s in all_fam[0:1] for s in freq_col_sample_names]

	print("family wise bit")
	print(familywise_bit_mask)

	# skip this iteration if none of infants are here 
	if not any(familywise_bit_mask[0:1]): #seed1 is this, seed2 is any(familywise_bit_mask[0])
		print("none of the infant samples are here")
		continue
	else:
		print("present")

	family_samples_present = [(s) for s in freq_col_sample_names if s in all_fam]
	print("family present")
	print(family_samples_present)

	if len(family_samples_present)  < 2:
		print("too little")
		continue

	# IMPOSE FILTER BASED ON READ DEPTH, snp_freq_matrix, snp_depth_matrix
	snp_depth_matrix_family =snp_depth_matrix[:,family_bit_mask]
	snp_freq_matrix_family = snp_freq_matrix[:,family_bit_mask]

	snp_freq_matrix_infant = snp_freq_matrix[:,infant_bit_mask]
	snp_depth_matrix_infant = snp_depth_matrix[:,infant_bit_mask]

	min_depth_bit_mask = snp_depth_matrix_family >= min_reads
	num_cols = snp_depth_matrix_family.shape[1]



	# if seed!= 6:
	# 	feature_filter = min_depth_bit_mask.sum(axis=1) >= num_cols 
	# else:
	feature_filter = np.where( (min_depth_bit_mask.sum(axis=1) >= num_cols) &( snp_freq_matrix_infant.sum(axis=1) > 0), True, False)
	print("LENGTH row names")
	print(len(freq_row_snp_names))
	freq_row_snp_names_family = np.array(freq_row_snp_names)[feature_filter]
	print("shape before filter")
	print(snp_freq_matrix_family.shape)
	print(snp_freq_matrix_family[0:10])
	print(snp_depth_matrix_family[0:10])
	snp_freq_matrix_family_filtered = snp_freq_matrix_family[feature_filter,:]
	snp_depth_matrix_family_filtered = snp_depth_matrix_family[feature_filter,:]



	# only taking last infant sample (so 12momnths)
	if study != "Li":
		last_infant_index = snp_freq_matrix_infant.shape[1]
	else:
		if seed == 1 or seed == 2:
			last_infant_index = 1 #Day 0 for odd seed, Day 2 for even seed
		if seed == 3 or seed == 4:
			last_infant_index = 2 # Day 2 for odd seed, Day 14 for even seed
		if seed == 5 or seed == 6:
			last_infant_index = 3 # Day 14 for odd seed, Day 42 for even seed
		if seed == 7 or seed == 8:
			last_infant_index = 4 # Day 42 for odd seed, Day 84 for even seed
		if seed == 9 :
			last_infant_index = 5 # Day 84 for odd seed

	snp_freq_matrix_infant_filtered = snp_freq_matrix_infant[feature_filter,(last_infant_index-1): last_infant_index]
	snp_depth_matrix_infant_filtered = snp_depth_matrix_infant[feature_filter,(last_infant_index-1): last_infant_index]

	print("snp_freq_matrix_family_filtered.shape")

	print(snp_freq_matrix_family_filtered.shape)
	print("len(freq_row_snp_names_family")
	print(len(freq_row_snp_names_family))


	print("check mother shape")
	snp_freq_matrix_filtered = snp_freq_matrix[feature_filter,:]
	snp_depth_matrix_filtered = snp_depth_matrix[feature_filter,:]
	snp_freq_matrix_mothers_filtered = snp_freq_matrix_filtered[:,mother_bit_mask]
	snp_depth_matrix_mothers_filtered = snp_depth_matrix_filtered[:,mother_bit_mask]
	print(snp_freq_matrix_mothers_filtered.shape)

	if custom_snv:

		snp_freq_matrix_family_filtered_private = snp_freq_matrix_family_filtered
		snp_depth_matrix_family_filtered_private = snp_depth_matrix_family_filtered
		private_snp_ind = list(range(snp_depth_matrix_family_filtered_private.shape[0]))
		
	else:
		likelihoods = score_snps(snp_freq_matrix_infant_filtered,snp_depth_matrix_infant_filtered,snp_freq_matrix_mothers_filtered, snp_depth_matrix_mothers_filtered)
		print("length likelihoods")
		print(len(likelihoods))

		# if study == "Dummy":
		# 	# print("ifant mother freq")
		# 	# print(snp_freq_matrix_infant_filtered)
		# 	# print(snp_freq_matrix_mothers_filtered)
		# 	# print("infant mother depth")
		# 	# print(snp_depth_matrix_infant_filtered)
		# 	# print(snp_depth_matrix_mothers_filtered)
		# 	# print(likelihoods)
		

		likelihood_for_cutoff = [l for l in likelihoods if l != -1]
		

		cutoff_best = 2*np.std(likelihood_for_cutoff) + np.mean(likelihood_for_cutoff) # sede 6 siulation
		print(cutoff_best)
		private_snp_ind = [ n for n,i in enumerate(likelihoods) if ((i > cutoff_best and i > 1) or i == -1) ]

		print("SD+2:  length private SNPs")
		print(len(private_snp_ind))

		if str(seed) == '99':
			private_snp_ind = random.sample(list(range(len(likelihoods))), k=len(private_snp_ind))
			print("Random:  length private SNPs")
			print(len(private_snp_ind))
			print(private_snp_ind)


		select_private_snps_bit_mask = [snp in private_snp_ind for snp in range(snp_freq_matrix_mothers_filtered.shape[0])]
		all_fam_LR[fam] = pd.Series(likelihoods )

		
		freq_row_snp_names_family = list(compress(freq_row_snp_names_family, select_private_snps_bit_mask))
		snp_freq_matrix_family_filtered_private = snp_freq_matrix_family_filtered[select_private_snps_bit_mask,:]
		snp_depth_matrix_family_filtered_private = snp_depth_matrix_family_filtered[select_private_snps_bit_mask,:]
		

	


	

	

	snp_counts_matrix_family_filtered_private = np.round_(np.multiply(snp_freq_matrix_family_filtered_private ,snp_depth_matrix_family_filtered_private),decimals = 0)
	ref_counts_matrix_family_filtered_private = np.round_(np.multiply((1-snp_freq_matrix_family_filtered_private) ,snp_depth_matrix_family_filtered_private),decimals = 0)

		# set to null no longer needed variables
	snp_freq_matrix_family_filtered = None
	snp_depth_matrix_family_filtered = None


	print(snp_freq_matrix_family_filtered_private)
	#print("filtered private row sum")
	#print(filtered_snp_freq_family_private.sum(axis=1)[0:10])
	# numpy_vers = False # numpy version is wrong
	if snp_counts_matrix_family_filtered_private.shape[0] > 0:
		print("private_snp_ind")
		print(private_snp_ind)

		family_private_snps_freq = pd.DataFrame(columns = all_fam,index=private_snp_ind)
		family_private_snps_freq.fillna(0)

		family_private_snps_counts = pd.DataFrame(columns = all_fam,index=private_snp_ind)
		family_private_snps_counts.fillna(0)

		family_private_refsnps_counts = pd.DataFrame(columns = all_fam,index=private_snp_ind)
		family_private_refsnps_counts.fillna(0)

		for c in range(len(family_samples_present)):
			family_private_snps_freq[family_samples_present[c]] = snp_freq_matrix_family_filtered_private[:,c]
			family_private_snps_counts[family_samples_present[c]] = snp_counts_matrix_family_filtered_private[:,c]
			family_private_refsnps_counts[family_samples_present[c]] = ref_counts_matrix_family_filtered_private[:,c]

		snp_freq_matrix_family_filtered_private = None
		snp_counts_matrix_family_filtered_private = None
		ref_counts_matrix_family_filtered_private = None
		final_family_private_snps_freq = family_private_snps_freq.values
		family_private_snps_freq_binned =final_family_private_snps_freq > familythreshold

		final_family_private_snps_counts = family_private_snps_counts.values
		final_family_private_refsnps_counts = family_private_refsnps_counts.values

		final_family_private_allsnps_counts = np.concatenate([final_family_private_snps_counts , final_family_private_refsnps_counts])

		#family_private_snps.to_csv(snp_dir +file_string + "_freq.csv")#, family_private_snps)#, delimiter=",")
		#family_private_snps_binned.to_csv(snp_dir +file_string + "_binned_freq.csv")#, family_private_snps_binned, delimiter=",")

		if not os.path.exists(out_snp_dir + "PrivateSNPs"):
			os.makedirs(out_snp_dir  + "PrivateSNPs")
		#np.savetxt(snp_dir + "PrivateSNPs/" + file_string + "_freq.csv", final_family_private_snps_freq, delimiter=",")
		family_private_snps_freq_binned = np.concatenate((np.transpose(np.array([freq_row_snp_names_family])),family_private_snps_freq_binned),axis=1)
		freq_row_snp_names_family1 = ["Alt_" + i for i in freq_row_snp_names_family]
		freq_row_snp_names_family2 = ["Ref_" + i for i in freq_row_snp_names_family]
		freq_row_snp_names_family = freq_row_snp_names_family1 + freq_row_snp_names_family2
		final_family_private_allsnps_counts = np.concatenate((np.transpose(np.array([freq_row_snp_names_family])),final_family_private_allsnps_counts),axis=1)

		print(final_family_private_allsnps_counts)
		np.savetxt(out_snp_dir  + "PrivateSNPs/" +file_string + "_binned_freq.csv", family_private_snps_freq_binned, delimiter=",",fmt='%s')
		np.savetxt(out_snp_dir  + "PrivateSNPs/" +file_string + "_counts.csv", final_family_private_allsnps_counts, delimiter=",",fmt='%s')

		print("--- %s seconds ---" % (time.time() - start_time))
	
	# END


if get_summary_stats:
	if not os.path.exists(out_snp_dir + "PrivateSNPs"):
		os.makedirs(out_snp_dir  + "PrivateSNPs")
	all_fam_LRs = pd.concat(all_fam_LR,axis=1)
	all_fam_LRs.to_csv(out_snp_dir + "PrivateSNPs/" +all_file_string + "_LRs.csv", sep=",")

	all_fam_presents = pd.concat(all_fam_present,axis=1)
	all_fam_presents.to_csv(out_snp_dir + "PrivateSNPs/" +all_file_string + "_fam_presence.csv",  sep=",")

