# author:    simon.wengert@helmholtz-muenchen.de
# date:      2024_04_23

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
    read_tsv("2023_05_17_gt_gPC_identified_Eur_donors.tsv") %>%
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
### ADD. CUSTOM FILTER TO REMOVE STRAND BIAS NOT ACCOUNTED BY MTDNA SERVER   ###
################################################################################

gt_rna_var_df <- gt_rna_var_df %>%
  # init strand bias metric, delta_frac, for highlighting this below
  mutate(delta_frac = abs(CoverageREV - CoverageFWD)/Coverage) %>%
  mutate(strand_bias = case_when(delta_frac > f_ops$delta_frac ~ TRUE, TRUE ~ FALSE))

#+ clean correlation plot of MAD outliers heteroplasmy calls vs couting --------
# Q: if we correlate the heteroplasmy levels, is there a clear pattern 
# for the potentially spurious ones (sum_major_minor != 1.00)?
# set MAD threshold
n_mad <- 10
# build plotting data frame
gt_rna_var_pre_strand_filter_df <- gt_rna_var_df %>%
  group_by(biospecimen_repository_sample_id,Pos) %>%
  mutate(sum_minor_level = sum(MinorLevel),
         # don't need to do sum for MajorLevel since these rows are just kept as
         # multiplets by mtDNA server
         sum_major_level = MajorLevel,
         sum_major_minor = sum_minor_level + sum_major_level,
         sum_major_minor_bigger_one = sum_major_minor > 1.00) %>%
  ungroup() %>%
  mutate(MAD = mad(sum_major_minor)) %>%
  mutate(mad_outlier = case_when(sum_major_minor >= 1.00 - MAD * n_mad & 
                                   sum_major_minor <= 1.00 + MAD * n_mad ~ FALSE,
                                 TRUE ~ TRUE))

# make plot pre strand filter
p_het_levels_mutserve_vs_counting_pre_filter <- gt_rna_var_pre_strand_filter_df %>%
  ggplot(aes(sum_heteroplasmic_level,sum_het_counts, colour = mad_outlier)) + 
  geom_point(size = 2, alpha = 0.6, stroke = 0.8, pch = 21) +
  geom_abline() +
  xlim(0, 1.00) +
  ylim(0, 1.00) +
  scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D"),   
                     guide = guide_legend(override.aes = list(shape = 15, alpha = 1,size = 6))) +
  p_ops$my_theme +
  facet_wrap(~molecular_process)
# make plot post strand filter
p_het_levels_mutserve_vs_counting_post_filter <- gt_rna_var_pre_strand_filter_df %>%
  filter(strand_bias == FALSE) %>%
  ggplot(aes(sum_heteroplasmic_level,sum_het_counts, colour = mad_outlier)) + 
  geom_point(size = 2, alpha = 0.6, stroke = 0.8, pch = 21) +
  geom_abline() +
  xlim(0, 1.00) +
  ylim(0, 1.00) +
  scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D"),   
                     guide = guide_legend(override.aes = list(shape = 15, alpha = 1,size = 6))) +
  p_ops$my_theme +
  facet_wrap(~molecular_process)
# assemple plot 
p_panel_strand_bias <-
  cowplot::plot_grid(p_het_levels_mutserve_vs_counting_pre_filter,
                     p_het_levels_mutserve_vs_counting_post_filter,
                     nrow = 2)

# save plot
ggsave(p_panel_strand_bias , 
       filename = '2024_04_28_p_het_levels_mutserve_vs_counting_custom_strand_bias_filter.png',
       width =  12,
       height = 12,
       dpi = "screen")   
ggsave(p_panel_strand_bias, 
       filename = '2024_04_28_p_het_levels_mutserve_vs_counting_custom_strand_bias_filter.pdf',
       width =  12,
       height = 12,
       dpi = 300)  

#+ get R2 values describing the filter ---------------
# fitting linear models to be able to get R2 values
lm_df <- gt_rna_var_pre_strand_filter_df
lm_dna_prior <- lm(sum_heteroplasmic_level ~ sum_het_counts, data = lm_df[which(lm_df$molecular_process == "dna_mutation") ,])
lm_dna_post <-  lm(sum_heteroplasmic_level ~ sum_het_counts, data = lm_df[which(lm_df$molecular_process == "dna_mutation" & lm_df$strand_bias == FALSE), ])
lm_rna_prior <- lm(sum_heteroplasmic_level ~ sum_het_counts, data = lm_df[which(lm_df$molecular_process == "rna_modification"),])
lm_rna_post <-  lm(sum_heteroplasmic_level ~ sum_het_counts, data = lm_df[which(lm_df$molecular_process == "rna_modification" & lm_df$strand_bias == FALSE), ])
# build data frame to plot/show those values
rsqr_df <- 
  tibble(rsqr_dna_prior = summary(lm_dna_prior)$r.squared,
         rsqr_dna_post = summary(lm_dna_post)$r.squared,
         rsqr_rna_prior = summary(lm_rna_prior)$r.squared, 
         rsqr_rna_post = summary(lm_rna_post)$r.squared) 
# print to console
rsqr_df

# A tibble: 1 × 4
# rsqr_dna_prior rsqr_dna_post rsqr_rna_prior rsqr_rna_post
#  <dbl>         <dbl>          <dbl>         <dbl>
#  0.945         0.977          0.920         0.919

#+ apply custom strand bias filer ----------------------------------------------
gt_rna_var_df <- gt_rna_var_df %>% filter(strand_bias == FALSE)
filter_numbers_df <- bind_rows(filter_numbers_df,
                               tibble(step = "strand_bias_custom",
                                      n_heteroplasmies = nrow(gt_rna_var_df),
                                      n_positions = length(unique(gt_rna_var_df$Pos)),
                                      n_samples = length(unique(gt_rna_var_df$biospecimen_repository_sample_id)),
                                      n_donors = length(unique(gt_rna_var_df$SUBJID))))


#+ reassign heteroplasmy levels using allelic counts for RNA modifications -----
gt_rna_var_df <- gt_rna_var_df %>%
  # keep a record of the original mtDNAserver heteroplasmy calls
  mutate(het_mutserve = heteroplasmic_level,
         sum_het_mutserve = sum_heteroplasmic_level) %>%
  mutate(heteroplasmic_level = as.double(heteroplasmic_level)) %>%
  mutate(heteroplasmic_level = case_when(molecular_process == "rna_modification" ~ het_counts,
                                         TRUE ~ heteroplasmic_level),
         sum_heteroplasmic_level =  case_when(molecular_process == "rna_modification" ~ sum_het_counts,
                                              TRUE ~ sum_heteroplasmic_level)) %>%
  dplyr::select(-strand_bias)

#+ make plot to show effects of reassignemnt -----------------------------------
p_replace_heteroplasmy_calls_with_allelic_counts_at_rna_mods <- gt_rna_var_df %>%
  ggplot(aes(sum_heteroplasmic_level,sum_het_counts)) + 
  geom_point(size = 2, alpha = 0.6, stroke = 0.8, pch = 21) +
  geom_abline() +
  xlim(0, 1.00) +
  ylim(0, 1.00) +
  p_ops$my_theme +
  facet_wrap(~molecular_process)

# save plot to disc
ggsave(p_replace_heteroplasmy_calls_with_allelic_counts_at_rna_mods, 
       filename = '2024_11_13_p_replace_heteroplasmy_calls_with_allelic_counts_at_rna_mods.pdf',
       width =  12,
       height = 6,
       dpi = 300)    

################################################################################
###       save filter-QCed heteroplasmies                                    ###
################################################################################

saveRDS(gt_rna_var_df,"2024_04_28_gt_rna_var_passing_het_filters_df.rds")
filter_numbers_df$filter_type <- "heteroplasmy_filter"
write_tsv(filter_numbers_df,"heteroplasmy_filter_numbers.tsv")


