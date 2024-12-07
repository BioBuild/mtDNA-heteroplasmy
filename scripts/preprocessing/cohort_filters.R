# author:    simon.wengert@helmholtz-muenchen.de
# date:      2024_04_28

library(tidyverse)
source("~/git/mtDNA_variants/scripts/utils/global_settings.R")
gt_rna_var_df <- readRDS("2024_04_28_gt_rna_var_n_het_min_positions_df.rds")
# load filter_number table to keep track of how much is remaining after each step
filter_numbers_df <- read_tsv("heteroplasmy_filter_numbers.tsv")

#+ apply filters ---------------------------------------------------------------
n_donors_per_tissue_df <- gt_rna_var_df %>%distinct(tissue,SUBJID) %>% count(tissue) %>% select(tissue, n_donors_per_tissue = n)
gt_rna_var_df <- 
    left_join(gt_rna_var_df,
              n_donors_per_tissue_df) %>%
    filter(n_donors_per_tissue >=  f_ops$tissue_n_donors_min)
# update filter numbers
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "tissue_n_donors_min",
                                      filter_type = "cohort_filter",
                                      n_heteroplasmies = nrow(gt_rna_var_df),
                                      n_positions = length(unique(gt_rna_var_df$Pos)),
                                      n_samples = length(unique(gt_rna_var_df$biospecimen_repository_sample_id)),
                                      n_donors = length(unique(gt_rna_var_df$SUBJID))))
# show filter numbers
filter_numbers_df

#+ save file to disc -----------------------------------------------------------
saveRDS(gt_rna_var_df,"2024_04_28_gt_rna_var_cohort_filters.rds")
# write out filter_numbers for commonality signal filter step
write_tsv(filter_numbers_df,"heteroplasmy_filter_numbers.tsv")
