# author:   simon.wengert@helmholtz-munich.de
# purpose:  load heteroplasmy genotype files, apply exact same filters as used 
#           in association testing, subset to RNA modifications and write out 
#           unique tissue-sample-pos lookup table for subsequent use in trans-
#           cript processing pipeline (based on pysam). The aim is to have the
#           read level summaries available for the same set of samples that went
#           into the association testing.

source("~/scripts/utils/global_settings.R")
data_path <- "/mtDNA_variants/"
model_input_path <- paste0(data_path,"metadata/genotype_files/split_by_tissue/")
tissues <- list.dirs(model_input_path, full.names = FALSE)[-1]
# load annotation of RNA modification sites 
mt_tRNA_modifications_df <- read_tsv("~/git/mtDNA_variants/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))

#+ load raw heteroplasmy estimates and allelic counts --------------------------
filename <- "2024_10_01_gt_rna_raw_position_lookup_table_rna_modifications.rds"
if(!file.exists(filename)){
  gt_rna_raw_df <- readRDS("2024_04_28_gt_raw_rna_heteroplasmy_estimates_allelic_counts_het_pos.rds")
  ## get sum coverage per Position
  gt_rna_raw_df <- gt_rna_raw_df %>%
    group_by(biospecimen_repository_sample_id,Pos) %>%
    mutate(sum_coverage = sum(Coverage)) %>%
    ungroup()
  ## subset to rna modifications
  gt_rna_raw_df <- gt_rna_raw_df %>% filter(Pos %in% m1A_G_methylations$genomic_pos)
  ## keep columns of interst
  gt_rna_raw_df <- gt_rna_raw_df %>% 
    distinct(Run,biospecimen_repository_sample_id,SUBJID,tissue,Pos,homoplasmic_base,sum_heteroplasmic_level,sum_coverage)
  ## save to disc
  saveRDS(gt_rna_raw_df,filename)
} else {
  gt_rna_raw_lookup_df <- readRDS(filename)
}

#+ loop through tissues and apply testing filters to each SNP ------------------
t_list <- list()
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("\n#tissues is ",t," out of #",length(tissues),"\n"))
  model_input_data <- read_tsv(paste0(model_input_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)
  
  # only consider mtRNA modification sites 
  model_input_data <- model_input_data[which(model_input_data$Pos %in% m1A_G_methylations$genomic_pos), ] 
  
  ################################################################################
  ###                           filters                                        ###
  ################################################################################

  #+ drop NAs for now doing complete case analysis -------------------------------
  model_input_data <- model_input_data %>% drop_na()

  #+ filter out positions 0 variance ---------------------------------------------
  model_input_data  <- 
    model_input_data %>%
    filter(var_per_pos > f_ops$site_variance_larger_than)
  
  #+ filter for minimum number of donors -----------------------------------------
  min_donors_pos <- 
    model_input_data  %>%
    group_by(Pos) %>%
    count(Pos) %>% 
    filter(n >= f_ops$site_n_donors_min) %>%
    ungroup() %>%  
    pull(Pos)
  
  # keep tissue table for selecting bam files later
  t_df <- model_input_data %>% 
    distinct(SUBJID,tissue,Pos,signal_type,sum_heteroplasmic_level,Coverage)
    
  # store in list for later use
  t_list[[t]] <- t_df
  
}
# collapse list to data frame
all_tissue_sample_pos_df <- bind_rows(t_list)
 
#+ add sample-tissue-pos level info including homoplasmic base -----------------
all_tissue_sample_pos_df <- 
  left_join(all_tissue_sample_pos_df,
            gt_rna_raw_lookup_df %>%
                distinct(Run,SUBJID,biospecimen_repository_sample_id,
                         tissue,Pos,homoplasmic_base),
            by = c("SUBJID","tissue","Pos"),relationship = "many-to-many") %>%
    distinct()

#+ save to disc ----------------------------------------------------------------
write_tsv(all_tissue_sample_pos_df,"2024_10_01_gt_sample_list_at_rna_modification_sites.tsv")
