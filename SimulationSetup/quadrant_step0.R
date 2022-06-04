# goal:
# (1 Dropout)
# Determine dropout rates in existing taxa (species_list)
# use beta draw
# (2 Complexity)
# already determined

# look back at empirical code
# get info for existing mother infant groupings in Seed 7 simulations
dropout_type = "high"
seed = 21
new_seed = 23

require(dplyr)
## EMPIRICAL TRANSMISSION CALCULATIONS - a mean - each species needs it's own list of families to test (test in get sigs)
data_dir = "/Users/leahbriscoe/Documents/FEASTX/"
ref_dir = paste0(data_dir,"BackhedFiles/")
new_dir = paste0(data_dir,"SimulationDFiles/")

# existing species list
species_list = read.csv(paste0(data_dir,
                               "code/pipeline/SimulationB_reference_info/species_list.txt"),sep="\t")

backhed_paper_data = read.csv(paste0(ref_dir,"metadata_merge.csv"),sep=",",stringsAsFactors = FALSE)
new_paper_data = read.csv(paste0(new_dir,"metadata_merge_seed", seed,".csv"),sep=",",stringsAsFactors = FALSE)
all_mom_acc = unlist(backhed_paper_data %>% filter(cohort == "M") %>% select(run_accession))
species_data = read.csv(paste0(ref_dir,"count_reads.txt"),sep="\t",stringsAsFactors = FALSE,row.names=1)
colnames(species_data)
species_data_rel = sweep(species_data, MARGIN = 2,STATS = colSums(species_data), "/")
species_data = species_data[,all_mom_acc]


species_data_rel = species_data_rel[,all_mom_acc]

##  Convert to function


# # gets probability based on multiple mother-birth pairs
# transmission_count = c()
# family_ids = unlist((backhed_paper_data %>% filter(cohort == "B") %>% select(STUDY.ID))$STUDY.ID)
# for(i in 1:length(family_ids)){
#   fam_data=backhed_paper_data %>% filter(STUDY.ID == family_ids[i]) 
#   infant_acc = unlist(fam_data %>% filter(cohort == "B") %>% select(run_accession))
#   mom_acc = unlist(fam_data %>% filter(cohort == "M") %>% select(run_accession))
#   transmission = species_data[,mom_acc] > 0 & species_data[,infant_acc] > 0
#   
#   colnames(species_data)
#   if(i == 1){
#     transmission_count = transmission
#   }else{
#     transmission_count = transmission_count + transmission
#   }
#   
# }
# dim(transmission_count)
# species_trans_prob = transmission_count/length(family_ids)
# names(species_trans_prob) = row.names(species_data)
# hist(species_trans_prob,breaks=100,
#      main="Over all 98 Backhed families, probability > 0 in B \n given > 0 in mom ")
#   
# hist(species_trans_prob[species_list[,1]],breaks=100,
#      main="[Just the 78 species with SNVs] \n Over all 98 Backhed families, probability > 0 in B \n given > 0 in mom ")


# Calculate high dropout transmission matrix of mothers
# Calculate low dropout transmission matrix of mothers
# lower probaility 0
low_dropout  = species_data_rel

for(row in 1:nrow(species_data_rel)){
  for(col in 1:ncol(species_data_rel)){
    b = 0.1 # low  dropout
    mean = species_data_rel[row,col]
    a = mean /(1-mean)
    
    low_dropout[row,col] = rbeta(n=1,shape1 = mean,shape2 = b)
    
  
  }
}

hist(unlist(rowMeans(low_dropout)),breaks=100,
     main="Low dropout: Mean transmission \n rate from beta a = relab(mom) ")

hist(unlist(rowMeans(low_dropout[species_list[,1],])),breaks=100,
     main="Low dropout:ç \n Mean transmission rate from beta a = relab(mom)")


high_dropout  = species_data_rel

for(row in 1:nrow(species_data_rel)){
  for(col in 1:ncol(species_data_rel)){
    b = 0.1 # low  dropout
    mean = species_data_rel[row,col] * 0.1
    a = mean /(1-mean)
    
    high_dropout[row,col] = rbeta(n=1,shape1 = mean,shape2 = b)
    
    
  }
}


hist(unlist(rowMeans(high_dropout)),breaks=100,
     main="High dropout: Mean transmission \n rate from beta a = relab(mom) ")

hist(unlist(rowMeans(high_dropout[species_list[,1],])),breaks=100,
     main="High dropout:\n Mean transmission rate from beta a = relab(mom)")



hist(rbeta(n = 1000, shape1 = 0, shape2 = 0.2),breaks=10)
hist(rbeta(n = 1000, shape1 = 0.005, shape2 = 0.1),breaks=100)


# then do binomial draw from each mother's beta distribution
## Make a high dropout dictionary for each study id. Label is complex or simple, 
## Make a low dropout dictionary for each study id. Label is complex or simple, 


#cor(species_trans_prob1[intersect_names],species_trans_prob2[intersect_names])


## ADD COMPLEXITY COLUMN
num_nonzero_sources = c()
complexity = c()
studies = unique(new_paper_data$study_id)
for(study_i in studies ){
  
  fam_data=new_paper_data %>% filter(STUDY.ID == study_i) 
  moms_nonzero = unlist(fam_data %>% filter(grepl("M",cohort), True_prop > 0) %>% select(run_accession))
  num_nonzero_sources = c( num_nonzero_sources,length(moms_nonzero))
  if(length(moms_nonzero) > 4){
    complexity  = c(complexity , "Complex")
    
  }else{
    complexity  = c(complexity , "Simple")
  }
  
}
complexity_meta = data.frame(STUDY.ID = studies, num_nonzero_sources,
                             complexity)

new_paper_data = merge(new_paper_data,complexity_meta,by = "STUDY.ID")
new_paper_data$dropout = dropout_type
write.table(new_paper_data ,paste0(new_dir,"metadata_merge_seed", new_seed,".csv"),
            row.names = TRUE,quote = FALSE,sep=",")



### RUN ONE BINOMIAL PER mom based on beta distirbutions
#input: dorpoutmatric
#output:species_presence_absence
if(dropout_type == "low"){
  dropout_matrix = low_dropout
}else{
  dropout_matrix = high_dropout
}

studies = unique(new_paper_data$study_id)
species_presence_absence = matrix(nrow = length(studies),ncol =  nrow(dropout_matrix))
row.names(species_presence_absence) = studies
colnames(species_presence_absence)  = row.names(dropout_matrix)
for(study_i in studies ){
  # study_i = studies[1]
  # species_i = "Alistipes_putredinis_61533"
  print(paste0("study_i ",study_i))
  
  fam_data=new_paper_data %>% filter(STUDY.ID == study_i) 
  infant_acc = unlist(fam_data %>% filter(cohort == "B") %>% select(run_accession))
  source_acc = unlist(fam_data %>% filter(True_prop > 0) %>% select(run_accession))
  # non zero mom only
  for(species_i in row.names(dropout_matrix)){
    
    transmission_vector = dropout_matrix[species_i,source_acc]
    yes_no_transmission_per_source = sapply(transmission_vector,  function(x){rbinom(1,1,prob = x)})
    if(any(yes_no_transmission_per_source > 0)){
      species_presence_absence[study_i,species_i] = 1
    }else{
      species_presence_absence[study_i,species_i] = 0
    }
    
  }
  
}

# transmission_vector = dropout_matrix["Alistipes_putredinis_61533",]
# sapply(transmission_vector,  function(x){rbinom(1,1,prob = x)})
# #plot number of positive present species per transmission
# species_presence_absence[,"Alistipes_putredinis_61533"]
hist(rowSums(species_presence_absence),breaks = 10,main = "All species number of \n considered species across families")

hist(rowSums(species_presence_absence[,species_list[,1]]),breaks = 10,main = "[Just the 78 species with SNVs] number of \n considered species across families")

write.table(species_presence_absence ,paste0(new_dir,"species_presence_absence_seed",new_seed,".csv"),
            row.names = TRUE,quote = FALSE,sep=",")

range(rowSums(species_presence_absence[,species_list[,1]]))
# from now on need to make sure species and SNvs considered are equivalent in terms of which species are considered

