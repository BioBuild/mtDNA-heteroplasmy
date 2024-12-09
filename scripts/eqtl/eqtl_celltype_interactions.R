# author:   simon.wengert@helmholtz-munich.de
library("tidyverse")
library("kableExtra")

genotype_path <- "/genotype_files/split_by_tissue/"
gene_expression_path <- "/gene_expression/"
save_path <- "/cis_eQTL/"
## read in xCell scores
cell_types_df <- read_tsv("2024_08_14_celltype_proportions_all_tissues.tsv")
cell_types_df <- cell_types_df %>% select(used_for_testing = enriched,everything())
cell_types_df <- cell_types_df %>% filter(used_for_testing == TRUE)
cell_types_df$proportion <- cell_types_df$int_xCell_score

# test in significant eQTLs  
gt_cis_eqtl_res_df <- read_tsv("/cis_eQTL/2024_07_23_gt_cis_eqtl_peer_factors_regressed_out_res_df.tsv")
eQTL_sig_df <- gt_cis_eqtl_res_df %>% filter(BB_E_bonferroni <= 0.05) %>% distinct()

## read in heteroplasmy data
gt_rna_var_df <- readRDS("2024_04_28_gt_rna_var_annotated.rds")

#+ perform interaction analysis per tissue -------------------------------------
# pull out tissues for iterating
tissues <- unique(cell_types_df$tissue)
t_res_list <- list()
for(t in seq_along(tissues)){   
  # build model input data per tissue
  t_rna_var_df <- gt_rna_var_df %>% filter(tissue == !! tissue) %>% distinct(biospecimen_repository_sample_id,SUBJID)
  t_cell_types_df <- cell_types_df[which(cell_types_df$tissue == tissue),]
  t_cell_types_df <- left_join(t_cell_types_df,t_rna_var_df) %>% distinct() 
  
  ## load gene expression files 
  t_ge_file <- read_tsv(paste0(gene_expression_path,tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"), show_col_types = FALSE)
  t_ge_file <- reshape2::melt(data.table::setDT(t_ge_file), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  t_ge_file <- left_join(t_ge_file,t_rna_var_df) %>% distinct() 
  
  ## load heteroplasmy genotypes
  t_genotypes_df <- read_tsv(paste0(genotype_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE) %>% 
    dplyr::select(snp_id = feature_id,everything()) %>%
    drop_na(sum_heteroplasmic_level)

  ## merge all 3 files from above into one model input data frame
  model_input_df <- left_join(t_genotypes_df,t_ge_file %>% select(-SUBJID), by = "biospecimen_repository_sample_id", relationship = "many-to-many") %>% distinct()
  model_input_df <- left_join(model_input_df,t_cell_types_df %>% select(-c("SUBJID","tissue","tissue_category")), by = "biospecimen_repository_sample_id", relationship = "many-to-many") %>% distinct()
    
  # proceed analysis only uisng complete cases
  n_samples_total <- length(unique(model_input_df$SUBJID))
  model_input_df <- model_input_df %>% drop_na(cell_type)
  n_samples_celltypes_available <- length(unique(model_input_df$SUBJID))
  model_input_df$n_samples_celltypes_available <- n_samples_celltypes_available
  model_input_df$frac_celltypes_available <- n_samples_celltypes_available/n_samples_total
  
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
    
    # apply test per cell_type
    cell_types <- as.character(unique(p_df$cell_type))
    c_res_list <- list()
    for(c in seq_along(cell_types)){ 
      # subset per cell_type
      cell_type <- cell_types[c]
      c_df <- p_df[which(p_df$cell_type == cell_type), ] 
      # only run test if variance is > 0 (i.e. if there is all 0 % for one c)
      # let's require 5 donors per cell type having more than 0 to be considered 
      if(var(c_df[ ,"proportion"]) > 0.00 && nrow(c_df[which(c_df$proportion == 0), ]) <= 5){
        
        ## interaction test inspired by xenophons deconvolution QTL
        ### define model formulas
        single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
        if(unique(p_df$tissue) %in% single_sex_tissues){   
          lm_m0_formula <- paste0("sum_heteroplasmic_level ~ logTPM + proportion + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
          lm_m1_formula <- paste0("sum_heteroplasmic_level ~ logTPM + proportion + logTPM : proportion + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        } else {
          lm_m0_formula <- paste0("sum_heteroplasmic_level ~ logTPM + proportion + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
          lm_m1_formula <- paste0("sum_heteroplasmic_level ~ logTPM + proportion + logTPM : proportion + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        }
        ### run linear models
        #### null model using additive effects
        lm_m0 <- lm(formula(lm_m0_formula), data = c_df)
        #### alternative model including cell type interaction term
        lm_m1 <- lm(formula(lm_m1_formula), data = c_df)
        res_m1 <- summary(lm_m1)
        #### comparing the result of the 2 models to see if adding an interaction term 
        #### adds more explainability to the the model -- using a likelihood ratio test
        LR_test <- lmtest::lrtest(lm_m1,lm_m0) 
        ### return result
        c_res_df <- data.frame(tissue = tissue,
                               pos_feature_id = pos_feature_id,
                               cell_type = cell_type,                 
                               n_samples = length(unique(t_genotypes_df$SUBJID)),
                               n_samples_cell_type_available = unique(c_df$n_samples_celltypes_available),
                               mean_cell_type_proportion = mean(c_df$proportion),
                               frac_celltypes_available = unique(c_df$frac_celltypes_available),
                               pass_testing_filters = TRUE,
                               lm_m1_interaction_estimate = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:proportion")), "Estimate"],
                               lm_m1_interaction_std_error = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:proportion")), "Std. Error"],
                               lm_m1_interaction_pval = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("logTPM:proportion")), "Pr(>|t|)"],
                               LR_pval = LR_test[2,'Pr(>Chisq)'])
        
      } else {
      # return NAs for those cell types not present or tested
      c_res_df <- data.frame(tissue = tissue,
                             pos_feature_id = pos_feature_id,
                             cell_type = cell_type,                 
                             n_samples = length(unique(t_genotypes_df$SUBJID)),
                             n_samples_cell_type_available = unique(c_df$n_samples_celltypes_available),
                             mean_cell_type_proportion = NA,
                             frac_celltypes_available = unique(c_df$frac_celltypes_available),
                             pass_testing_filters = FALSE,
                             lm_m1_interaction_estimate = NA,
                             lm_m1_interaction_std_error = NA,
                             lm_m1_interaction_pval = NA,
                             LR_pval = NA)  
        
      }
      # save to list
      c_res_list[[c]] <- c_res_df
      
    }
    # save to list
    p_res_df <- bind_rows(c_res_list)
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
res_df <- res_df %>%
  filter(pass_testing_filters == TRUE) %>%
  mutate(LM_A_bonferroni =  p.adjust(LR_pval, method = 'bonferroni')) %>%
  dplyr::select(-pass_testing_filters)


#+ plot significant results to console -----------------------------------------
res_df %>%
  filter(LM_A_bonferroni <= 0.05) %>%
  kbl() %>%
  kable_paper(bootstrap_options = "striped", full_width = F)


#+ save to disc ----------------------------------------------------------------
write_tsv(res_df,paste0(save_path,"2024_08_14_LM_A_cis_eQTL_all_tissues_cell_type_deconvolution.tsv"))

