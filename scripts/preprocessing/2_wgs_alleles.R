# author:    simon.wengert@helmholtz-muenchen.de
# date:      2023_05_16
library(tidyverse)
library(data.table)
source("~/scripts/utils/global_settings.R")
################################################################################
### read in WGS bam file derived mtDNA server calls.                         ###
################################################################################

#+ read in mtDNA_server raw outputs and apply filters and QC -------------------
wgs_raw_filename_v9 <- paste0("/wgs/2023_07_06_gt_raw_wgs_samples_df.rds")
if(!file.exists(wgs_raw_filename_v9)){

# wgs samples v9: 
wgs_v9_var_dir <- "/wgs/mtvar/"
wgs_v9_filenames <- list.files(path = wgs_v9_var_dir , 
                                pattern = "*.txt",
                                full.names = FALSE) 
wgs_v9_filenames <- wgs_v9_filenames[grep("_raw.txt",wgs_v9_filenames)]
# length(wgs_v9_filenames)
# [1] 897
# length(unique(wgs_v9_filenames))
# [1] 897
#
# read in mtDNA_server variant calls & put into tidy df format
ldf <- list()
for(i in seq_along(wgs_v9_filenames)){  
    # load data
    ldf[[i]] <- fread(paste0(wgs_v9_var_dir,wgs_v9_filenames[i]))
    # rescue biospecimen_repository_sample_id from input filename
    ldf[[i]]$Run <- strsplit(wgs_v9_filenames[i], c("_raw.txt"))       
    ldf[[i]]$file_count <- paste0("file_nr_",i)
}
wgs_samples_v9 <- 
  rbindlist(ldf) %>% 
  mutate(biospecimen_repository_sample_id = Run) %>%
  dplyr::select(file_name = SAMPLE,biospecimen_repository_sample_id,everything()) 

# apply filters and define homoplasmic base - AKA inherited allele
# stash <- wgs_samples_v9
wgs_samples_v9 <- 
    wgs_samples_v9 %>%
        dplyr::select(top_fwd = 'TOP-FWD',
                      cov_fwd = 'COV-FWD',
                      top_rev = 'TOP-REV',
                      cov_rev = 'COV-REV',
                      everything()) %>%
        # maybe this one is too stringent
        filter(cov_fwd > as.integer(f_ops$min_cov_fwd), 
               cov_rev > as.integer(f_ops$min_cov_rev)) %>%
        # define homoplasmic_base
        mutate(homoplasmic_base = case_when(top_fwd == top_rev ~ top_fwd,
                                            cov_fwd > cov_rev ~ top_fwd,
                                            cov_fwd < cov_rev ~ top_rev)) %>%
           dplyr::select(SUBJID,POS,REF,top_fwd,top_rev,cov_fwd,cov_rev,homoplasmic_base,
                      tissue,COHORT,SEX,AGE,RACE,ETHNCTY,HGHT,WGHT,BMI) 
# dim(wgs_samples_v9)
# [1] 14854155       17
# length(unique(wgs_samples_v9$SUBJID))
# [1] 868

saveRDS(wgs_samples_v9,wgs_raw_filename_v9)
} else {
    gt_wgs_raw_var_df <- readRDS(wgs_raw_filename_v9)
}
# dim(gt_wgs_raw_var_df)
# [1] 14861436       17

################################################################################
###  add WGS derived inherited_allele to rnaseq heteroplasmies calls         ###
################################################################################

gt_rna_var_df <- 
  left_join(gt_rna_var_df,
            gt_wgs_raw_var_df %>%
              dplyr::select(SUBJID,Pos = 'POS',inherited_allele = homoplasmic_base),
              by = c('SUBJID','Pos')) %>% 
  distinct()
# dim(gt_rna_var_df)
# [1] 1454773      32

# keep WGS base at Pos where available - fill in REF base otherwhise
gt_rna_var_df <-
	gt_rna_var_df %>%
        mutate(covered_by_wgs  = case_when(is.na(inherited_allele) ~ FALSE,
										 TRUE ~ TRUE)) %>%
		    mutate(inherited_allele = case_when(is.na(inherited_allele) ~ Ref,
			 							 TRUE ~ inherited_allele))

# remove rna heteroplasmies without WGS coverage ------------------------------
gt_rna_var_df <- gt_rna_var_df %>% filter(covered_by_wgs == TRUE)
table(gt_rna_var_df$covered_by_wgs)
#  TRUE 
# 1696595 
length(unique(gt_rna_var_df$SUBJID))
# [1] 840
# this means we have RNA-seq derived mt-vars for 840 donors

################################################################################
### update homoplasmic allele using inherited allele obtained in WGS         ###
################################################################################

gt_rna_var_df$heteroplasmic_base <- ifelse(gt_rna_var_df$inherited_allele == gt_rna_var_df$MajorBase, 
                                           gt_rna_var_df$MinorBase, 'NA')
gt_rna_var_df$heteroplasmic_base <- ifelse(gt_rna_var_df$inherited_allele == gt_rna_var_df$MinorBase,
                                           gt_rna_var_df$MajorBase, gt_rna_var_df$heteroplasmic_base)

gt_rna_var_df$heteroplasmic_base <- ifelse(gt_rna_var_df$inherited_allele != gt_rna_var_df$MinorBase &
                                           gt_rna_var_df$inherited_allele != gt_rna_var_df$MajorBase,
                                           gt_rna_var_df$MinorBase, gt_rna_var_df$heteroplasmic_base) 
#
# define heteroplasmic_level
gt_rna_var_df$heteroplasmic_level <- ifelse(gt_rna_var_df$inherited_allele == gt_rna_var_df$MajorBase,
                                            gt_rna_var_df$MinorLevel, 'NA')
gt_rna_var_df$heteroplasmic_level <- ifelse(gt_rna_var_df$inherited_allele == gt_rna_var_df$MinorBase,
                                            gt_rna_var_df$MajorLevel, gt_rna_var_df$heteroplasmic_level)
gt_rna_var_df$heteroplasmic_level <- ifelse(gt_rna_var_df$inherited_allele != gt_rna_var_df$MinorBase &
                                            gt_rna_var_df$inherited_allele != gt_rna_var_df$MajorBase,
                                            gt_rna_var_df$MinorLevel, gt_rna_var_df$heteroplasmic_level) 

# and assign a value to heteroplasmic_level. 
gt_rna_var_df$heteroplasmic_level <- ifelse(gt_rna_var_df$MinorBase == "-" &
                                            gt_rna_var_df$inherited_allele == gt_rna_var_df$Variant,  
                                            1.000 - gt_rna_var_df$VariantLevel, 
                                            gt_rna_var_df$heteroplasmic_level)
gt_rna_var_df$heteroplasmic_level <- ifelse(gt_rna_var_df$MinorBase == "-" &
                                            gt_rna_var_df$inherited_allele != gt_rna_var_df$Variant,  
                                            gt_rna_var_df$VariantLevel, 
                                            gt_rna_var_df$heteroplasmic_level)        

# assigning homoplasmic allele
gt_rna_var_df$homoplasmic_base <- gt_rna_var_df$inherited_allele

saveRDS(gt_rna_var_df,"/2023_07_06_gt_rna_before_flipping_df.rds")
