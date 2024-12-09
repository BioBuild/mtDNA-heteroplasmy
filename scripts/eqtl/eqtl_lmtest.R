# author:  simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library(lmtest)
metadata_dir <- "/genotype_files/split_by_tissue/"
gene_expression_path <- "/gene_expression/"
results_dir <- "/cis_eQTL/"
tissues <- list.dirs(gene_expression_path, full.names = FALSE)[-1]

#+ loop through tissues and fit lm model per tissue to find candidate SNPs -----
LM_cis_eQTL_A_start_time <- Sys.time()
t_res_list <- list()
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("\n#tissues is ",t," out of #",length(tissues),"\n"))
  # load genotypes
  # model_input_data <- read_tsv(paste0(model_input_path,tissue,"/2023_06_06_gt_",tissue,"_model_donor_age_relationship_input_data.tsv"))
  t_genotypes <- 
    read_tsv(paste0(metadata_dir,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE) %>% 
    select(snp_id = feature_id,everything())
  positions <- unique(t_genotypes$Pos)
  # filter out rows with 0 variance 
  # apply **testing stage filters** to genotypes
  input_data <- t_genotypes
  # 1) drop samples with missing values
  input_data <- input_data %>% drop_na()
  # 2) filter out positions 0 variance
  input_data <- input_data %>% filter(var_per_pos > f_ops$site_variance_larger_than)
  # 3) filter for min number of donors 
  min_donors_pos <- 
    input_data  %>%
    group_by(Pos) %>%
    count(Pos) %>% 
    filter(n >= f_ops$site_n_donors_min) %>%
    ungroup() %>%  
    pull(Pos)
  positions <- positions[which(positions %in% min_donors_pos)]
  #
  # load gene expression
  t_ge_file <- read_tsv(paste0(gene_expression_path,tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"), show_col_types = FALSE)
  # make long to wide
  t_ge_file <- melt(setDT(t_ge_file), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  # join gene expression files to genotypes
  input_data <- left_join(input_data,t_ge_file, by = "SUBJID")
  # loop through each mtDNA position and test for
  # pairwise association with genes
  features <- unique(input_data$feature_id)
  # fitting the model
  pos_res_list <- list()
  for(p in seq_along(positions)){
    position <- positions[p]
    pos_input_data <- 
      input_data %>%
      filter(Pos == position)
    feature_res_list <- list()
    for(f in seq_along(features)){
      feature <- features[f]
      feature_input_data <- 
        pos_input_data %>%
        filter(feature_id == feature)
      # fit cis-eQTL
      lm_m0 <- lm(sum_heteroplasmic_level ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = feature_input_data)
      lm_m1 <- lm(sum_heteroplasmic_level ~ logTPM + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = feature_input_data)
      res_m1 <- summary(lm_m1)
      LR_test <- lrtest(lm_m1,lm_m0) # go on here
      tmp_df <- data.frame(tissue = tissue,
                           Pos = position,
                           feature_id = feature,
                           n_samples = length(unique(feature_input_data$SUBJID)),
                           estimate =  lm_m1$coefficients["logTPM"],
                           lm_std_error = res_m1$coefficients["logTPM","Std. Error"],
                           LR_pval =  LR_test[2,  'Pr(>Chisq)'])
      feature_res_list[[f]] <- tmp_df
      cat(paste0("\n#positions ",p,"/",length(positions)," #features ",f,"/",length(features),"\n"))
    }
    pos_res_list[[p]] <- bind_rows(feature_res_list)
  }
  # collapse list
  t_res_df <- bind_rows(pos_res_list)
  # return result
  t_res_list[[t]] <- t_res_df
}
# collapse list
lm_eqtl_analytical_df <- bind_rows(t_res_list)

#+ record total run time -------------------------------------------------------
LM_cis_eQTL_A_end_time <- Sys.time()
lm_eqtl_analytical_df$run_time_hours <- abs(round(difftime(LM_cis_eQTL_A_start_time,LM_cis_eQTL_A_end_time, units = "hours"), digits = 5))

#+ subset to polyA genes ---------------------------------------------------
lm_eqtl_analytical_df <-
  left_join(lm_eqtl_analytical_df,
             g_ops$mt_regions %>% dplyr::select(feature_id = ID ,poly_A),
             by = "feature_id")
lm_eqtl_analytical_df <- lm_eqtl_analytical_df %>% filter(poly_A == TRUE)
# multiple testing correction study wide
lm_eqtl_analytical_df$analytical_FDR <- p.adjust(lm_eqtl_analytical_df$LR_pval, method = 'fdr')
lm_eqtl_analytical_df$analytical_bonferroni = p.adjust(lm_eqtl_analytical_df$LR_pval, method = 'bonferroni')

#+ save to disc ----------------------------------------------------------------
write.table(lm_eqtl_analytical_df,paste0(results_dir,"2024_07_23_gt_cis_eQTL_LM_A_all_tissues_poly_A_peer_factors_regressed_out.tsv"))

