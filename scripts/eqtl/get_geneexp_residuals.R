# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
covariate_path <- "/covariates/" ## peer factors 
gene_expression_path <- "/gene_expression/"
tissues <- dir(gene_expression_path)

## remove kidney cortex from analysis as too low sample size <70
tissues <- tissues[-which(tissues == "kidney_cortex")]

## regress out tissue-specific peer factors 
# loop through tissues
for(t in seq_along(tissues)){
  # set tissue names
  t_small_caps <- tissues[t]
  t_large_caps <- as.character(gt_lookup_small_vs_large_caps_tissues_names_df[t,"tissues_xeno"])
  # tell us which tissue is being processed
  cat(paste0("tissue #",t,"/",length(tissues),"\n"))
  # load peer factor file
  t_peer_factors_df <- read_tsv(paste0(covariate_path,t_large_caps,"_allcovs.txt"))
  t_peer_factors_df$SUBJID <- t_peer_factors_df$sample_id
  t_peer_factors_df <- t_peer_factors_df[,c("SUBJID",grep("PEER",colnames(t_peer_factors_df),value = TRUE))]
  # load gene expression file
  t_ge_logTMP_polyA_wide_df <- read_tsv(paste0(gene_expression_path,t_small_caps,"_mtgenes_poly_A.txt"))
  t_ge_logTMP_polyA_long_df <- melt(setDT(t_ge_logTMP_polyA_wide_df), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  # merge gene expression file with PEER factors
  t_ge_peer_df <- left_join(t_ge_logTMP_polyA_long_df,t_peer_factors_df,by = "SUBJID")
  # regressing out PEER factors from gene expression files
  ## pull out peer factors
  peer_factors <- grep("PEER",colnames(t_ge_peer_df),value = TRUE)
  cat(paste0("# peer factors is ",length(peer_factors),"\n"))
  ## specify model formula
  model_formula <- paste0("logTPM ~ ",paste0(peer_factors,collapse = " + "))
  ## run linear regression to regress out peer factors from gene expression
  res_lm <- lm(formula(paste(model_formula, collapse = " ")), data = t_ge_peer_df)
  ## get the residuals of that regression and keep them as corected logTPMs
  t_ge_peer_df$logTPM_corrected <- res_lm$residuals
  # make wide version of th 
  t_ge_peer_factors_regressed_out_logTPM_polyA_wide_df <- t_ge_peer_df %>% 
    dplyr::select(feature_id,SUBJID,logTPM_corrected) %>% 
    pivot_wider(names_from = SUBJID,values_from = logTPM_corrected)
  # save to corrected gene expression file to disc
  write_tsv(t_ge_peer_factors_regressed_out_logTPM_polyA_wide_df,
            paste0(gene_expression_path,t_small_caps,"_mtgenes_poly_A_peer_factors_regressed_out.txt"))
}

#+ qc reporting ----------------------------------------------------------------
write_tsv(qc_df,"2024_06_13_logTPM_regress_out_peer_factors_correlation_before_after.tsv")
