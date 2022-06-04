args = commandArgs(trailingOnly=TRUE)

# args =  c("snps_only",1, 10, "1.0", "Infant", 10, "Brooks" ,0 )
# local_bool = TRUE
# plot_local =TRUE
local_bool = FALSE # local_bool = TRUE
plot_local =FALSE # plot_local =TRUE


debug_mode = FALSE
filter_unique_SNV = TRUE
#
inclusion_criteria =args[1]#"combo"#"otus_only" #"combo"#snps_only"#"otus_only" # 
start=as.integer(args[2])
print(args)
seed=as.integer(args[3])
uniqueness=1
threshold=args[4]
sink = args[5]
min_reads = args[6]
study = args[7]
dropout = as.integer(args[8])
familythreshold=0.1
select_genus = FALSE
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
lapply(Packages, library, character.only = TRUE)

sink_alt = sink
if(sink == "4M"){
  sink_alt = "X4M"
}else if (sink == "12M"){
  sink_alt = "X12M"
}

#Rscript SourceTrackingScript.R otus_only 1 10 1.0 B 10 SimulationD -1 
require(data.table)
## END packages

if(local_bool){
  input_dir =  paste0('~/Documents/FEASTX/',study,'Files/')
  code_dir = '~/Documents/FEASTX/code/pipeline/'
  overall_dir = '~/Documents/FEASTX'
}else if(grepl("Suez",study)){
  input_dir =  paste0('/u/scratch/b/briscoel/',study,'Files/')
  code_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/code/pipeline/'
  overall_dir = '/u/scratch/b/briscoel/'
}else{
  input_dir =  paste0('/u/home/b/briscoel/project-ngarud/FEASTX/',study,'Files/')
  code_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/code/pipeline/'
  overall_dir = '/u/home/b/briscoel/project-ngarud/FEASTX/'
}
if(((seed >= 10) & (seed <= 13)) & study == "SimulationD"){
  seed_data = 7
}else if((seed %in% c(115,215,109,209)) & study == "Tara"){
  seed_data = (seed  %% 100)
}else if((seed==22 | seed == 23) & study == "SimulationD"){
  seed_data = 21
}else{
  seed_data = seed
}
equal_sampling_1d_4d = FALSE

#library("FEAST-master", lib.loc = paste0(overall_dir))
### PREPARING COUNT DATA
counts_ab = fread(paste0(input_dir,"count_reads.txt"),sep="\t",stringsAsFactors = FALSE)
save_rows = counts_ab$species_id
counts_ab = as.matrix(counts_ab[,2:ncol(counts_ab)])
row.names(counts_ab) = save_rows
print("original counts ab dim")
print(dim(counts_ab))  
# SNV_species = read.csv(paste0(input_dir,"species_list_seed",seed,".txt"),stringsAsFactors = FALSE,header=FALSE)
# 
# if(grepl("otu",inclusion_criteria)){
#   
#   intersection = intersect(SNV_species[,1],row.names(counts_ab))
#   print("length interesction")
#   print(length(intersection))
#   counts_ab = counts_ab[intersection,,drop=FALSE]
#   counts_ab_orig = counts_ab
#   
# }
counts_ab_orig = counts_ab

### PREPARING METADATA
fam_reference = read.csv(paste0(input_dir,"Pipe_ChosenFamily_seed",seed,".csv"),stringsAsFactors = FALSE)
if(study == "SimulationD"){
  meta_merge = read.csv(paste0(input_dir,"metadata_merge_seed",seed,".csv"),stringsAsFactors = FALSE)
  
}
fam_ids = sort(unique(fam_reference$family_id))
file_string = paste0("_seed_",seed,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                     familythreshold,"_minreads_",min_reads)

inputfile_string = paste0("_seed_",seed_data,"_uniqueness_",uniqueness,"_threshold_",threshold,"_familythresh_",
                     familythreshold,"_minreads_",min_reads)


if (dropout < 0){
  transmission_0_1 = read.csv(paste0(input_dir,"species_presence_absence_seed",seed,".csv"),
                                  row.names = 1,stringsAsFactors = FALSE)

  
}

set.seed(0)
ALLtimes= list()

for(fam_id_i in c(start:length(fam_ids))){
  # fam_id_i = 48
  ##################  SELECT SAMPLES: based on arrangement ###########
  print("fam number")
  print(fam_id_i)
  fam_id = fam_ids[fam_id_i]
  print("Fam id")
  print(fam_id)
  fam_ref = fam_reference %>% filter(family_id == fam_id)
  select_samples = fam_ref[2:length(fam_ref)]
  if(study == "Metasub"){
    names_select_samples = names(select_samples)[select_samples != "haib17DB4959_H3MGVCCXY_SL259919"]
    select_samples = select_samples[select_samples != "haib17DB4959_H3MGVCCXY_SL259919"]
    names(select_samples) = names_select_samples
  }
  print("number sources")
  print(length(select_samples))
  
  
  
  ##################  END SELET       ################## 
  
  
  ################### FILTER OUT SPECIES ##########
  if (dropout < 0){
    non_zero_species = names(transmission_0_1[fam_id,transmission_0_1[fam_id,] != 0])
    print("non zero species")
    print(non_zero_species)
  }
  
  if(grepl("otu",inclusion_criteria) &  (dropout < 0)){
    intersect_transmission = intersect(row.names(counts_ab_orig),non_zero_species)
    print(counts_ab_orig[1:4,1:4])
    counts_ab = counts_ab_orig[intersect_transmission,,drop=FALSE]
    print("Dim(counts_ab")
    print(dim(counts_ab))
  }
  
    
  ################### LOAD: SNPs ##################
  if(inclusion_criteria == "snps_only"){
    path_check = paste0(input_dir,"CATTED_fam_",fam_id,inputfile_string,"_counts.csv.bz2")
    
    if(equal_sampling_1d_4d){
      prefix = "D"
    }else{
      prefix = "ED"
    }
    if(((seed == 12) & study == "SimulationD") | ((seed %in% c(109,115)) & study == "Tara")){
      path_check = paste0(input_dir,prefix,"1_CATTED_fam_",fam_id,inputfile_string,".csv.bz2")
    }
    if(((seed == 13) & study == "SimulationD") | ((seed %in% c(209,215)) & study == "Tara")){
      path_check = paste0(input_dir,prefix,"4_CATTED_fam_",fam_id,inputfile_string,".csv.bz2")
      
    }
    
    if(!file.exists(file.path(path_check))){
      print(path_check)
      print(paste0("Family file ", fam_id, " absent"))
      next
    }

    
    
    print(path_check)
    
    matrix_dat <- tryCatch(fread(path_check,stringsAsFactors = FALSE,header = FALSE,
                                 quote ="",na.strings=c("NA","NaN", " ","nan","")), error=function(e){
                                   print("not file")
                                   return(e)
                                 })
    if(inherits(matrix_dat, "error")){
      ALLtimes[[fam_id]] = NA
      next
    }else{
      matrix_dat =as.matrix(matrix_dat)
      print(dim(matrix_dat))
    }
    print(matrix_dat[1:4,1:4])
    keep_name =  matrix_dat[,1]
    
    mat_index = which(fam_ref %in% select_samples)
    matrix_dat = matrix_dat[,mat_index]
    matrix_dat = apply(matrix_dat, 2,as.integer)
    matrix_dat[is.na(matrix_dat)] = 0
    colnames(matrix_dat) = select_samples
    row.names(matrix_dat) = keep_name
    
    # subsample so sampel SNV number
    
    # if(((seed == 12 ) | (seed ==13)) & study == "SimulationD"){
    #   
    #   if(file1_length != file2_length){
    #     print('original number snvs')
    #     print(nrow(matrix_dat))
    #     
    #     # OUTDATED LINE
    #     #draw = sample(1:nrow(matrix_dat), size=min_snvs, replace=FALSE)
    #     #matrix_dat = matrix_dat[draw,]
    #     
    #     
    #     ## FIXING 
    #     
    #     just_species = sapply(row.names(matrix_dat),function(x){
    #       temp = gsub("Ref_","",gsub("Alt_","",x))
    #       return(temp)
    #     })
    #     print("length unique just species should be half")
    #     print(length(just_species))
    #     print(length(unique(just_species)))
    #     draw = sample(unique(just_species), size=(min_snvs/2), replace=FALSE)
    #     matrix_dat = matrix_dat[(just_species %in% draw),]
    #     
    #     
    #     print('new number snvs')
    #     print(nrow(matrix_dat))
    #   }
    #   
    # }
    # 
    if(dropout < 0){
      rnames = row.names(matrix_dat)
      rnames_spt = strsplit(rnames,"\\|")
      just_species = do.call(rbind,rnames_spt)[,1]
      just_species = sapply(just_species,function(x){
        temp = gsub("Ref_","",gsub("Alt_","",x))
        return(temp)
      })
      chosen_snps = sapply(just_species,function(x){
        if(x %in% non_zero_species){
          return(TRUE)
        }else{
          return(FALSE)
        }
      })
      print(sum(chosen_snps))
      print("dim before dropout")
      print(dim(matrix_dat))
      matrix_dat = matrix_dat[chosen_snps,]
      print("dim after dropout")
      print(dim(matrix_dat))
    }
  }
  
  if(inclusion_criteria == "otus_only"){
    print("Select_samples")
    print(select_samples)
    print("colnames")
    print(colnames(counts_ab))
    print(intersect(colnames(counts_ab),select_samples))
    intersectselect_samples = intersect(colnames(counts_ab),as.character(select_samples))
    matrix_dat =  counts_ab[,as.character(select_samples)]
  }
  
  #  MAKE METADATA FOR SOURCE TRACK
  col_id = rep(fam_id,length(select_samples))
  col_env = names(select_samples)
  col_sourcesink = rep("Source",length(select_samples))
  col_sourcesink[which(col_env== sink_alt)] = "Sink"
  
  print("colenv")
  print(col_env)
  print("col_siourcessink")
  print(col_sourcesink)
  print("colid")
  print(col_id)
  custom_metadata = data.frame('Env' = col_env, 'SourceSink' = col_sourcesink, 'id' = col_id)
  row.names(custom_metadata) = col_env
  colnames(matrix_dat) = col_env
  print(custom_metadata)
  t1 = Sys.time()  
  

 meta_dat = custom_metadata #LB 4/14/22
  
  ## LB meta_dat = custom_metadata[colSums(matrix_dat)!=0,] 
  print("Zero column check")
  #print(matrix_dat[1:20,])
  print(colSums(matrix_dat))
  if(sum(colSums(matrix_dat)!=0) == 1){
    print("skippy")
    next
  }
  ## LB matrix_dat = matrix_dat[,colSums(matrix_dat)!=0]
  #dim(matrix_dat)
  
  print("before feast input")
  print(dim(matrix_dat))
  if(dim(matrix_dat)[1] == 0){
    next
  }
  print('sum NA')
  print(sum(is.na(matrix_dat)))
  
  
  
  
  #matrix_dat = matrix_dat[rowSums(matrix_dat)!=0,]
  
  # only non zero samples
  if(study == "SimulationD"){
    meta_fam = meta_merge %>% filter(study_id == fam_id)
    print("meta fam")
    print(meta_fam)
    
  }
  print("checking first lines")
  
  if(debug_mode){
    overall_dir = '~/Documents/FEASTX'
    source(paste0(overall_dir,"FEAST/R/FEAST.R"))
    FEAST_output <- tryCatch(FEAST(C = t(matrix_dat), metadata =meta_dat, different_sources_flag = 0, dir_path =input_dir,
                                   outfile="demo",overall_dir =paste0(overall_dir,"FEAST/R/")), error=function(e){
                                     print(e)
                                     print("feast not work")
                                     return(e)})
  }else{
    #install.packages("devtools")
    #install.packages("rlang")
    #
    
    #dim(matrix_dat)
    library(FEAST)
    
    if(filter_unique_SNV){
      print("CHECKING FOR UNIQUE SNV")
      print(dim(matrix_dat))
      time1 = Sys.time()
      vector = c(1,rep(0,ncol(matrix_dat)-1))
      bool_matrix = t(t(matrix_dat) == vector)
      print("total length of time")
      print(Sys.time() - time1)
      unique_check = bool_matrix[,1]
      
      # get row names of "unique" to remove alt and ref
      unique_snv = row.names(matrix_dat)[unlist(unique_check)]
      unique_snv = gsub("Alt_","",gsub("Ref_","",unique_snv))
      all_snv = gsub("Alt_","",gsub("Ref_","",row.names(matrix_dat)))
      bool_vector = !(all_snv %in% unique_snv)
      matrix_dat = matrix_dat[bool_vector,]
      print(dim(matrix_dat))
    }
    
    
    #dim(meta_dat)
    # test = t(matrix_dat)
    # metadata =meta_dat
    # 
    colsum_sources = colSums(matrix_dat[,2:ncol(matrix_dat)])
    csum = colSums(matrix_dat)
    coverage_min =min(csum[csum!=0])
   
    if(sum(colsum_sources > 1) == 1){
      matrix_dat = rbind(matrix_dat,c(0,rep(1,(ncol(matrix_dat) -1))))
    }
    
   
    
    FEAST_output <- tryCatch(FEAST(C = t(matrix_dat), metadata =meta_dat, different_sources_flag = 0, dir_path =input_dir,
                                   outfile="demo",COVERAGE =coverage_min), error=function(e){
                                     print(e)
                                     print("feast not work")
                                     return(e)})

  }
  

  
  if(inherits(FEAST_output, "error")){
    #REAL WORK
    ALLtimes[[fam_id]] = NA
    #next
  }else{
    print("FEAST TEST")
    FEAST_output = FEAST_output[c(1:(nrow(meta_dat)))]
    orig_names_output = names(FEAST_output)
    FEAST_output = do.call(rbind,FEAST_output)
    print(Sys.time() - t1)
    all_times = data.frame(variable =  orig_names_output , value = FEAST_output[,1])
    print(all_times)
    
    ALLtimes[[fam_id]] = all_times
  }
  
  if (equal_sampling_1d_4d){
    saveRDS(ALLtimes,paste0(input_dir,"EqFEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".rds"))
    
  }else{
    saveRDS(ALLtimes,paste0(input_dir,"FEAST_box_input_", inclusion_criteria,"_start_", start,file_string, "_dropout_", dropout,"_sink_",sink,".rds"))
    
  }
  
  
  
}


