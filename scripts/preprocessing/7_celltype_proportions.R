# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library("tidyverse")
library("kableExtra")

#+ load and QC-check GTEx provided xCell scores used in sarahs paper -----------
# GTEx_Analysis_v8_xCell_scores_7_celltypes.txt has been downloaded from 
# https://gtexportal.org/home/downloads/adult-gtex/qtl on 2024-06-27
file_name <- "2024_07_26_GTEx_Analysis_v8_xCell_scores_7_celltypes_long_form_annotated_by_simon.txt"
if(!file.exists(file_name)){
  gt_xCell_scores_df <- read_tsv("~/metadata/annotations/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt", show_col_types = FALSE)
  # lookuptable for annotating tissues, SUBJIDs etc.v
  gt_v8_lookup_df <- read_tsv("~/metadata/annotations/gt_v8_lookup_df.tsv", show_col_types = FALSE)
  # wrangle and annotate
  gt_xCell_scores_df <- gt_xCell_scores_df %>%
    pivot_longer(-c("cell_type"),names_to = "biospecimen_repository_sample_id",values_to = "xCell_enrichment_score")
  gt_xCell_scores_df <- 
    left_join(gt_xCell_scores_df,
              gt_v8_lookup_df %>% 
                dplyr::select(biospecimen_repository_sample_id,tissue,tissue_category,SUBJID)) %>% 
    distinct()
  # annotate with large cap names for GTEx nomencalture
  gt_xCell_scores_df <- 
    left_join(gt_xCell_scores_df,
              gt_v8_lookup_df %>% 
                dplyr::select(biospecimen_repository_sample_id,tissue,tissue_category,SUBJID)) %>% 
    distinct()
  # gt_xCell_scores_df
  gt_xCell_scores_df <- 
    left_join(gt_xCell_scores_df,
            gt_lookup_small_vs_large_caps_tissues_names_df %>%
              transmute(tissue = tissues_simon, tissues_xeno))
  # annotate if present in our set of RNAseq heteroplasmy calls 
  gt_xCell_scores_df <- gt_xCell_scores_df %>%
    mutate(present_in_rna_var = case_when(biospecimen_repository_sample_id %in% unique(gt_rna_var_df$biospecimen_repository_sample_id) ~ TRUE,
                                          TRUE ~ FALSE))
    # save file to disc for Xenofon for applying transformations
  write_tsv(gt_xCell_scores_df,paste0("~/metadata/annotations/",file_name))
} else {
  gt_xCell_scores_df <- read_tsv(paste0("~/metadata/annotations/",file_name))
}

#+ xCell score threshold used by Kim-Hellmuth et al. 2020 ----------------------
xCell_thr <- 0.1 
# caclculate median xCell score 
gt_xCell_scores_df <- gt_xCell_scores_df %>%
  group_by(tissue,cell_type) %>%
    mutate(t_median_xCell_enrichment_score = median(xCell_enrichment_score)) %>%
  ungroup()
# annotate if passing xCell_thr 
gt_xCell_scores_df <- gt_xCell_scores_df %>%
  mutate(enriched = case_when(t_median_xCell_enrichment_score > 0.1 ~ TRUE,
                              TRUE ~ FALSE))

#+ apply inverse normal transformation per tissue and cell type ----------------
inverse_normal_transform <- function(x) {
  ranks <- rank(x, ties.method = "average")
  uniform_quantiles <- ranks / (length(ranks) + 1)
  normal_scores <- qnorm(uniform_quantiles)
  return(normal_scores)
}
if(F){
  int_function <- function(x){
    qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  }
}
# apply inverse normal transformation
gt_xCell_scores_df <- 
  gt_xCell_scores_df %>%
    group_by(cell_type,tissue) %>%
      mutate(int_xCell_score = inverse_normal_transform(xCell_enrichment_score)) %>%
      # mutate(int_xCell_score = int_function(xCell_enrichment_score)) %>%
    ungroup()

#+ save all tissues cell type proportion file to disc --------------------------
write_tsv(gt_xCell_scores_df, "~/metadata/annotations/2024_08_14_celltype_proportions_all_tissues.tsv")


#+ load heteroplasmy files and compare the number of donor ids ---------
gt_rna_var_df <- readRDS("~/2024_04_28_gt_rna_var_annotated.rds")
gt_rna_var_tissue_subject_df <- gt_rna_var_df %>%
  distinct(tissue,SUBJID,biospecimen_repository_sample_id) %>%
  filter(!tissue == "kidney_cortex")
# make n_donors file
n_donors_per_tissue_df <- gt_rna_var_tissue_subject_df %>% 
  group_by(tissue) %>% 
  summarise(n_donors_heteroplasy_genotypes = dplyr::n()) %>% 
  ungroup()
# check cell type file availability
cell_type_availability_df <- 
  left_join(gt_rna_var_tissue_subject_df,
            cell_types_df %>% 
              select(SUBJID,tissue,n_donors_cell_type_file) %>%
              mutate(cell_types_available = TRUE),
            by = c("SUBJID","tissue")) %>%
  mutate(cell_types_available = case_when(is.na(cell_types_available) ~ FALSE, 
                                          TRUE ~ cell_types_available)) %>%
  distinct()
# how many do we have cell type info available?
# table(cell_type_availability_df$cell_types_available)
# FALSE  TRUE 
# 3431  9101 
# check xCell score file numbers and include as well 
gt_xCell_scores_present_summary_df <- gt_xCell_scores_df %>%
  filter(biospecimen_repository_sample_id %in% unique(gt_rna_var_df$biospecimen_repository_sample_id)) %>%
  distinct(biospecimen_repository_sample_id,tissue) %>%
  count(tissue) %>%
  transmute(tissue,n_xCell_score_file = n)
# are these a subset?
length(which(unique(gt_rna_var_df$biospecimen_repository_sample_id) %in% unique(gt_xCell_scores_df$biospecimen_repository_sample_id) == TRUE))
# yes, these are a complete subset
# [1] 12532
