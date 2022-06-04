# SNV-FEAST


### Accessions for shotgun metagenomic data downloads

Backhed et al. 2015 (Mother infant): PRJEB6456	
Brooks et al. 2017 (NICU): PRJEB323631

Sunagawa et al. 2015 (Tara Oceans): PRJEB402


### Needed MIDAS files for signature SNV selection for FEAST Source Tracking


For each species with sufficient coverage:

* ref_freq.txt 
* snp_depth.txt

### Signature SNV selection

run Step_1.py to make appriprate input files needed for signature SNV selection and FEAST source tracking later on. 

run Step_2.py to run signature SNV selection. By default the signature SNV selection is done in 200 kbp windows

run Steps 3 through 5 to merge and zip files appropriately. 


### FEAST source tracking

input: metadata file produced by Step 1 above.
input: private SNVs file produced by the Signature SNV selection

** For simulation** Additional input needed





