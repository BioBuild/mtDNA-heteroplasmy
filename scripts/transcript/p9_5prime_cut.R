# authors:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library(lmtest)
library(gridExtra)
library("kableExtra")
# filter param
n_reads_modified_min <- 10
n_samples_per_pos_min <- 60

#+ source aux files ------------------------------------------------------------
# gtex sample lookup file
gt_v8_lookup_df <- read_tsv("gt_v8_lookup_df.tsv", show_col_types = FALSE)
# donor gPCs
gt_gPCs <- read_tsv("gt_v8_gPCs_europ.tsv")
# RNA modification sites 
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv", show_col_types = FALSE)
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))
p9_positions <- m1A_G_methylations %>% filter(transcript_pos == 9) %>% pull(genomic_pos)
# load RNA modification sites previously tested by Ali et al.
pos_gene_df <- read_tsv("2024_11_11_positions_tested_in_Ali_et_al_2020.tsv")

#+ load data -------------------------------------------------------------------
all_tissues_summary_df <-
  read_tsv("2024-10-08_all_tissues_sample_summary_from_read_level_analysis_filtered.tsv",
           show_col_types = FALSE) %>%
  filter(Pos %in% p9_positions)

if(F){
#+ restrict analysis to significant results from previous replication ----------
all_tissues_summary_df <- all_tissues_summary_df %>%
  left_join(.,replication_results_df %>% select(tissue,Pos,sig_read_counts)) %>%
  filter(sig_read_counts == TRUE)
}

#+ restrict analysis to RNA modfications considered by Ali et al ---------------
# reduction from  62345 to 40,732 (~ 65 %)
all_tissues_summary_df <- all_tissues_summary_df %>%
  left_join(.,pos_gene_df %>% mutate(keep = TRUE) %>% select(Pos,keep)) %>%
  filter(keep == TRUE) %>%

#+ annotate tissue summaries with donor covariates -----------------------------
all_tissues_model_input_df <- all_tissues_summary_df %>%
  mutate(n_reads_not_modified = n_reads_passing - n_reads_modified) %>%
  select(biospecimen_repository_sample_id,tissue,Pos,feature_strand,
         heteroplasmic_level_mtDNA_server,read_counts_heteroplasmic_level,n_reads_passing,
         n_reads_not_modified,n_reads_modified,n_is_modified_T_cut_five_prime_T,
         # lets keep in some additional rows but only for contrasting the 
         # trends in the visualization of the results. Not for modeling!
         n_is_modified_T_cut_five_prime_F,n_is_modified_F_cut_five_prime_T,
         n_is_modified_F_cut_five_prime_F) %>%
  left_join(.,gt_v8_lookup_df %>% distinct(biospecimen_repository_sample_id,SUBJID,SEX,AGE)) %>%
  left_join(.,gt_gPCs) %>%
  select(biospecimen_repository_sample_id,SUBJID,everything()) %>%
  distinct()
  select(-keep)

if(F){
# only test in tissues > 60 samples 
all_tissues_model_input_df <- 
  all_tissues_model_input_df %>% 
  group_by(tissue,Pos) %>%
   mutate(n_samples_per_pos = dplyr::n()) %>%
  ungroup() %>%
  filter(n_samples_per_pos >= n_samples_per_pos_min)
}

#+ perform position level inference for methylation and cutting ----------------
# do this for all samples looking at individual cut sites
tissues <- unique(all_tissues_model_input_df$tissue)
t_list <- list()
for(t in seq_along(tissues)){
  cat(paste0("tissue #",t,"/",length(tissues)," --------------------------\n"))
  # subset to tissue of interest and 
  tissue <- tissues[t]
  # keep only reads that are modified and not overlapping cut site
  t_df <- all_tissues_model_input_df %>% filter(tissue == !!tissue) 
  #loop through positions now
  positions <- unique(as.character(t_df$Pos))
  p_list <- list()
  for(p in seq_along(positions)){
    # subset to postion 
    position <- positions[p]
    p_df <- t_df %>% filter(Pos == !! position)
    # keep track of number of samples
    n_samples <- length(unique(p_df$biospecimen_repository_sample_id))
    
    single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
    if(unique(p_df$tissue) %in% single_sex_tissues){   
      null_model_formula <- c("cbind(n_is_modified_T_cut_five_prime_T,n_is_modified_T_cut_five_prime_F) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(n_is_modified_T_cut_five_prime_T,n_is_modified_T_cut_five_prime_F) ~ read_counts_heteroplasmic_level + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    } else {
      null_model_formula <- c("cbind(n_is_modified_T_cut_five_prime_T,n_is_modified_T_cut_five_prime_F) ~ SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(n_is_modified_T_cut_five_prime_T,n_is_modified_T_cut_five_prime_F) ~ read_counts_heteroplasmic_level + SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }
    
    # Initialize warning tracking
    null_warning <- FALSE
    alt_warning <- FALSE
    
    # Fit binomial null model with warning tracking
    null_model <- tryCatch({
      glm(as.formula(null_model_formula), family = "binomial", data = p_df)
    }, warning = function(w) {
      null_warning <<- TRUE
      glm(as.formula(null_model_formula), family = "binomial", data = p_df)
    })
    
    # Fit binomial alternative model with warning tracking
    alt_model <- tryCatch({
      glm(as.formula(alt_model_formula), family = "binomial", data = p_df)
    }, warning = function(w) {
      alt_warning <<- TRUE
      glm(as.formula(alt_model_formula), family = "binomial", data = p_df)
    })
    
    # summarise results for extracting info
    alt_model_summary <- summary(alt_model)
    
    # likelihood ratio test for assessing significance
    LR_test <- lrtest(null_model,alt_model)
    
    # assemble result
    tmp_df <- 
      data.frame(tissue = tissue,
                 Pos = position,
                 feature_strand = unique(p_df$feature_strand),
                 n_samples = n_samples,
                 mean_n_reads_modified = unique(p_df$mean_n_reads_modified),
                 # model specs and results
                 model_type = "binomial",
                 model_formula = alt_model_formula,
                 estimate = alt_model_summary$coefficients["read_counts_heteroplasmic_level","Estimate"],
                 std_error = alt_model_summary$coefficients["read_counts_heteroplasmic_level", "Std. Error"],
                 p_value = LR_test$`Pr(>Chisq)`[2],
                 null_warning = null_warning,
                 alt_warning = alt_warning)

    #  return result
    p_list[[p]] <- tmp_df
  }
  # collapse list and return results
  t_res_df <- bind_rows(p_list)
  t_list[[t]] <- t_res_df
}
# collapse list
res_df <- bind_rows(t_list)
# do multiple testing correction
## study wide bonferroni
res_df$bonferroni <- p.adjust(res_df$p_value,method = "bonferroni")
## study wide BH
res_df$BH <- p.adjust(res_df$p_value,method = "BH")
## within tissue (as this seems closest to what Ali et al. had done.)
res_df <- res_df %>%
  group_by(tissue) %>%
  mutate(t_bonferroni = p.adjust(p_value, method = "bonferroni"),
         t_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() 


#+ save results to disc --------------------------------------------------------
write_tsv(res_df, "2024_11_11_res_position_level_regression_five_prime.tsv")

