# author:    simon.wengert@helmholtz-muenchen.de
# date:      2023_07_06
library(tidytext)
source("~/scripts/utils/global_settings.R")
all_filters_df <- read_csv("~/metadata/utils/all_filters.csv")
gt_rna_var_df <- readRDS("2024_02_14_gt_rna_var_pre_filtering_het_df.rds")

gt_rna_var_df  <- 
  gt_rna_var_df %>%
    distinct() %>%
    group_by(biospecimen_repository_sample_id,Pos) %>%
      mutate(sum_heteroplasmic_level = sum(as.double(heteroplasmic_level))) %>%
    ungroup()

gt_high_heteroplasmies <- 
    gt_rna_var_df %>% 
        filter(sum_heteroplasmic_level > 0.5) %>%
        dplyr::select(SUBJID,tissue,Pos,Ref,Variant,homoplasmic_base,sum_heteroplasmic_level,covered_by_wgs,biospecimen_repository_sample_id)
# nrow(gt_high_heteroplasmies)
# [1] 107632
# how many of the heteroplasmies are bigger than 50 % in general
# nrow(gt_high_heteroplasmies)/nrow(gt_rna_var_df)
# [1] 0.06344001 --> around 0.6 %

################################################################################
### filter for EUR         donors                                            ###
################################################################################

gt_eurpean_donors <- 
    read_tsv("~/metadata/annotations/2023_05_17_gt_gPC_identified_Eur_donors.tsv") %>%
    pull(SUBJID)
gt_rna_var_df <- gt_rna_var_df %>% filter(SUBJID %in% gt_eurpean_donors)

#+ order columns of gt_rna_var_df nicely ---------------------------------------
colnames_in_nice_order <-
  c("SUBJID","Pos","Ref","homoplasmic_base","heteroplasmic_base",
    "heteroplasmic_level","sum_heteroplasmic_level","Variant","VariantLevel",
    "MajorBase","MajorLevel","MinorBase","MinorLevel","Coverage","tissue",
    "tissue_category","AGE","SEX","PMI","RACE","ETHNCTY","HGHT","WGHT","BMI")
gt_rna_var_df <- 
  gt_rna_var_df %>%
    dplyr::select(colnames_in_nice_order,everything())

################################################################################
### identify mtRNA modifications                                             ###
################################################################################

mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))
gt_rna_var_df <- gt_rna_var_df %>% mutate(molecular_event = case_when(Pos %in% m1A_G_methylations$genomic_pos ~ "rna_modification",
                                                                      TRUE ~ "dna_mutation"))

################################################################################
###               REMOVE MTDNA SITES PRONE TO SEQ ARTEFACTS                  ###
################################################################################
#+ apply f_ops heteroplasmy filter thr as defined in utils/all_filters csv -----
# and exlude microsattelites flagged by Peter Campbell's lab
gt_rna_var_df <- gt_rna_var_df %>% filter(!Pos %in% f_ops$sites_seq_artefacts)
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "seq_artefacts",
                                      n_heteroplasmies = nrow(gt_rna_var_df)))

################################################################################
###                 ONLY SELECT HETEROPLASMIC VARIANTS                       ###
################################################################################
gt_rna_var_df <- gt_rna_var_df %>% filter(Type == f_ops$variant_type)
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "variant_type",
                                      n_heteroplasmies = nrow(gt_rna_var_df)))

################################################################################
###                       filter for strand bias                             ###
################################################################################
gt_rna_var_df <- gt_rna_var_df %>% filter(Filter == f_ops$mutserve_filter)
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "mutserve_filter_pass",
                                      n_heteroplasmies = nrow(gt_rna_var_df)))
# just for reference: this can have the following values (+ we loose a lot)
#     PASS        STRAND_BIAS 
#     298528      738532 

################################################################################
###                       FILTER FOR COVERAGE                                ###
################################################################################
gt_rna_var_df <- 
  gt_rna_var_df %>%
  filter(CoverageFWD > as.integer(f_ops$min_cov_fwd), 
         CoverageREV > as.integer(f_ops$min_cov_rev))
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "min_coverage",
                                      n_heteroplasmies = nrow(gt_rna_var_df)))

################################################################################
###       save filter-QCed heteroplasmies                                    ###
################################################################################

saveRDS(gt_rna_var_df,"2024_04_04_gt_rna_var_qced_pre_wmh_filtering_df.rds")


