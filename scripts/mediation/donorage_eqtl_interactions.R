# author:   simon.wengert@helmholtz-munich.de
library("tidyverse")
library("kableExtra")

gt_rna_var_df <- readRDS("2024_04_28_gt_rna_var_annotated.rds")
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))

save_path <- "/mediation_analysis/"
genotype_path <- "/genotype_files/split_by_tissue/"
gene_expression_path <- "/gene_expression/"

# let's only use significant eQTLs now 
gt_cis_eqtl_res_df <- read_tsv("/cis_eQTL/2024_07_23_gt_cis_eqtl_peer_factors_regressed_out_res_df.tsv")
eQTL_sig_df <- gt_cis_eqtl_res_df %>% filter(BB_E_bonferroni <= 0.05) %>% distinct()

#+ perform interaction analysis per tissue -------------------------------------
tissues <- unique(gt_rna_var_df$tissue)
t_res_list <- list()
for(t in seq_along(tissues)){ 
  # subset for tissue
  tissue <- tissues[t]
  cat(paste0("tissue #",t,"/",length(tissue),"\n"))
  
  # build model input data per tissue
  t_rna_var_df <- gt_rna_var_df %>% filter(tissue == !! tissue) %>% distinct(biospecimen_repository_sample_id,SUBJID)
  
  ## load gene expression files 
  t_ge_file <- read_tsv(paste0(gene_expression_path,tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"), show_col_types = FALSE)
  ## make long to wide
  t_ge_file <- reshape2::melt(data.table::setDT(t_ge_file), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  ## annotate with biospecimen sample ID
  t_ge_file <- left_join(t_ge_file,t_rna_var_df, by = c("SUBJID")) %>% distinct() 
  
  ## load heteroplasmy genotypes
  t_genotypes_df <- read_tsv(paste0(genotype_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE) %>% 
    dplyr::select(snp_id = feature_id,everything()) %>%
    drop_na(sum_heteroplasmic_level)
  ## merge all 2 files from above into one model input data frame
  model_input_df <- left_join(t_genotypes_df,t_ge_file %>% select(-SUBJID), by = "biospecimen_repository_sample_id", relationship = "many-to-many") %>% distinct()

  # subset to pos_feature for which we found significant eQTL
  t_pos_feature_df <- eQTL_sig_df[which(eQTL_sig_df$tissue == tissue & eQTL_sig_df$BB_E_bonferroni <= 0.05), c("Pos","feature_id")]
  t_pos_feature_df$pos_feature_id <- paste0(t_pos_feature_df$Pos,"_",t_pos_feature_df$feature_id) 
  pos_feature_ids <- unique(t_pos_feature_df$pos_feature_id)
  model_input_df$pos_feature_id <- paste0(model_input_df$Pos,"_",model_input_df$feature_id) 
  model_input_df <- model_input_df %>% filter(pos_feature_id %in% pos_feature_ids)
  
  # apply test below to each eQTL position feature pair
  p_res_list <- list()
  for(p in seq_along(pos_feature_ids)){
    # subset for pos_feature pairs
    pos_feature_id <- pos_feature_ids[p]
    p_df <- model_input_df[which(model_input_df$pos_feature_id == pos_feature_id), ] 
    
    ## interaction test inspired by xenophons deconvolution QTL
    ### define model formulas
    single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
    if(unique(p_df$tissue) %in% single_sex_tissues){   
      lm_m0_formula <- paste0("sum_heteroplasmic_level ~ logTPM + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      lm_m1_formula <- paste0("sum_heteroplasmic_level ~ logTPM + AGE + logTPM : AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    } else {
      lm_m0_formula <- paste0("sum_heteroplasmic_level ~ logTPM + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      lm_m1_formula <- paste0("sum_heteroplasmic_level ~ logTPM + AGE + logTPM : AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }
    ### run linear models
    #### null model using additive effects
    lm_m0 <- lm(formula(lm_m0_formula), data = p_df)
    #### alternative model including cell type interaction term
    lm_m1 <- lm(formula(lm_m1_formula), data = p_df)
    res_m1 <- summary(lm_m1)
    #### comparing the result of the 2 models to see if adding an interaction term 
    #### adds more explainability to the the model -- using a likelihood ratio test
    LR_test <- lmtest::lrtest(lm_m1,lm_m0) 
    ### return result
    p_res_df <- data.frame(tissue = tissue,
                           pos_feature_id = pos_feature_id,
                           n_samples = length(unique(p_df$SUBJID)),
                           lm_m1_interaction_estimate = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:AGE")), "Estimate"],
                           lm_m1_interaction_std_error = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:AGE")), "Std. Error"],
                           lm_m1_interaction_pval = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:AGE")), "Pr(>|t|)"],
                           LR_pval = LR_test[2,'Pr(>Chisq)'])
    
    # save to list
    p_res_list[[p]] <- p_res_df
  }
  # save to list
  t_res_df <- bind_rows(p_res_list)
  t_res_list[[t]] <- t_res_df
}

# bind rows to get final results data frame
res_df <- bind_rows(t_res_list)
res_df <- res_df %>% separate(pos_feature_id, c("Pos","feature_id"))

#+ do study-wide bonferroni correction -----------------------------------------
res_df <- res_df %>% mutate(LM_A_bonferroni =  p.adjust(LR_pval, method = 'bonferroni')) 


#+ plot significant results ---------------------------------------------------
res_df %>%
  filter(LM_A_bonferroni <= 0.05) %>%
  mutate(molecular_process = case_when(Pos %in% m1A_G_methylations$genomic_pos ~ "rna_modification",
                                       TRUE ~ "dna_mutation")) %>%
  kbl() %>%
  kable_paper(bootstrap_options = "striped", full_width = F)


#+ plot QQ for sanity check ----------------------------------------------------
gap::qqunif(res_df$LR_pval)


#+ save to disc ----------------------------------------------------------------
write_tsv(res_df,paste0(save_path,"2024_09_23_LM_A_donor_age_logTPM_interaction_for_cis_eQTL_sig_hits.tsv"))

