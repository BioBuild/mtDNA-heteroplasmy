library(tidyverse)
library(data.table)
source("~/scripts/utils/global_settings.R")
##############################################################################
###      load mtDNA_server rna seq variant calls                           ###
##############################################################################

mutserve_var_dir <- "/mtvar/"
variant_filenames <- list.files(path = mutserve_var_dir, 
                                pattern = "*.txt",
                                full.names = FALSE) 
variant_filenames <- variant_filenames[-grep("_raw.txt",variant_filenames)]

## quick check for duplicates
length(variant_filenames)
# [1] 12740
length(unique(variant_filenames))
# [1] 12740
# --> nice, no duplicates 

## read in all variant calls 
ldf <- list()
for(i in seq_along(variant_filenames)){  
    # load data
    ldf[[i]] <- fread(paste0(mutserve_var_dir,variant_filenames[i]))
    ldf[[i]]$Run <- strsplit(variant_filenames[i], c(".txt"))       
    ldf[[i]]$file_count <- paste0("file_nr_",i)
}

## combine data from all samples 
samples_v8 <- 
    rbindlist(ldf) %>% 
      mutate(biospecimen_repository_sample_id = Run) %>%
      dplyr::select(file_name = ID, biospecimen_repository_sample_id,everything()) 

## get metadata 
gt_v8_public_biospecimen_repository_sample_id <- 
  read_tsv("~/metadata/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") %>%
    filter(SMAFRZE == "RNASEQ") %>%
    dplyr::select(biospecimen_repository_sample_id = 'SAMPID', everything()) %>%
    pull(biospecimen_repository_sample_id)

# FILTER our RNA-seq derived heteroplasmies to GTEx v8 samples only:
gt_rna_var_df <- 
    samples_v8 %>%
        filter(biospecimen_repository_sample_id %in% gt_v8_public_biospecimen_repository_sample_id)

##############################################################################
###      annotate cohort metadata & reference trinucleotide profiles       ###
##############################################################################


#+ annotate gt_rna_var_df with trinucleotide profiles ------------------------
lookup_trinuc_context_df <- 
  read_tsv("~/metadata/annotations/mtDNA_ref_trinucleotide_profiles_lookup_df.tsv")
gt_rna_var_df <- 
  gt_rna_var_df %>%
    left_join(.,lookup_trinuc_context_df %>% dplyr::select(Pos, ref_trinuc_profile))
dim(gt_rna_var_df)
# [1] 1881694      18

##############################################################################
###                        checking for duplicates                         ###
##############################################################################


#+ safetycheck for duplicates ------------------------------------------------
# Q: are there duplicated filenames?
duplicates_df <- 
    gt_rna_var_df %>%
        distinct() %>% 
        distinct(SUBJID,tissue,file_name) %>%
        group_by(SUBJID,tissue) %>%
            summarise(count_duplicates = dplyr::n()) %>%
        ungroup() %>% 
        filter(count_duplicates > 1) %>%  
        arrange(desc(count_duplicates))
# dim(duplicates_df)
# [1] 0 3
# YES!!! PERFECT --> no duplicates in there

##############################################################################
###                        save to disk                                    ###
##############################################################################

# save to disc
  saveRDS(gt_rna_var_df,rna_raw_filename_v9)
} else {
  gt_rna_var_df <- readRDS(rna_raw_filename_v9)
}

