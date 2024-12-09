# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library("tidyverse")
library("kableExtra")

save_path <- "/donor_age/"
genotype_path <- "/genotype_files/split_by_tissue/"
cell_types_df <- read_tsv("2024_08_14_celltype_proportions_all_tissues.tsv", show_col_types = FALSE) %>%
  mutate(used_for_testing = enriched, proportion = int_xCell_score) %>%
  distinct(SUBJID,biospecimen_repository_sample_id,tissue,cell_type,proportion,xCell_enrichment_score,int_xCell_score,used_for_testing)

plot_pheno_res_df <- read_tsv("2024_07_23_gt_donor_age_all_tissues_annotated.tsv", show_col_types = FALSE)
gt_xCell_scores_df <- read_tsv("2024_07_26_GTEx_Analysis_v8_xCell_scores_7_celltypes_long_form_annotated_by_simon.txt", show_col_types = FALSE)

## 1) keep only tissue-cell type combination with median xCell >= 0.1
xCell_thr <- 0.1 
### caclculate median xCell score 
gt_xCell_scores_df <- gt_xCell_scores_df %>%
  group_by(tissue,cell_type) %>%
  mutate(t_median_xCell_enrichment_score = median(xCell_enrichment_score)) %>%
  ungroup()
### annotate if passing xCell_thr 
gt_xCell_scores_df <- gt_xCell_scores_df %>%
  mutate(enriched = case_when(t_median_xCell_enrichment_score > 0.1 ~ TRUE,
                              TRUE ~ FALSE))
### keep only those which are passing the threshold
gt_xCell_scores_df <- gt_xCell_scores_df %>% filter(enriched == TRUE)

## 2) apply inverse normal transformation
inverse_normal_transform <- function(x) {
  ranks <- rank(x, ties.method = "average")
  uniform_quantiles <- ranks / (length(ranks) + 1)
  normal_scores <- qnorm(uniform_quantiles)
  return(normal_scores)
}
gt_xCell_scores_df <- 
  gt_xCell_scores_df %>%
    group_by(cell_type,tissue) %>%
      mutate(int_xCell_score = inverse_normal_transform(xCell_enrichment_score)) %>%
    ungroup()

cell_types_df$proportion <- cell_types_df$int_xCell_score

#+ perform donor age celltype interaction analysis per tissue ------------------
tissues <- unique(cell_types_df$tissue)
data_to_test_df <- donor_age_sig_df
# loop through tissues
tissues <- tissues[tissues %in% unique(data_to_test_df$tissue)]
t_res_list <- list()
for(t in seq_along(tissues)){ 
  # subset for tissue
  tissue <- tissues[t]
  cat(paste0("tissue #",t,"/",length(tissues),"\n"))
  # build model input data per tissue
  ## load and format cell type table
  t_cell_types_df <- cell_types_df[which(cell_types_df$tissue == tissue) ,]
  ## load heteroplasmy genotypes
  t_genotypes_df <- read_tsv(paste0(genotype_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE) %>% 
    dplyr::select(snp_id = feature_id,everything()) %>%
    drop_na(sum_heteroplasmic_level)
  ## merge all 2 files from above into one model input data frame
  model_input_df <- 
    left_join(t_genotypes_df, 
              cell_types_df %>% 
                distinct(biospecimen_repository_sample_id,tissue,
                         cell_type,proportion),
                         by = c("biospecimen_repository_sample_id","tissue"),
                         relationship = "many-to-many")
  #-----------------------------------------------------------------------------
  ## merge all 2 files from above into one model input data frame
  #model_input_df <- left_join(t_genotypes_df,t_cell_types_df, by = c("SUBJID","tissue"),relationship = "many-to-many") # HERE CHECK MATCH MAKING AS YOU MAY CAN DO THAT WITH SAMPLE IDS TOO.
  
  # proceed analysis only using complete cases
  n_samples_total <- length(unique(model_input_df$SUBJID))
  model_input_df <- model_input_df %>% drop_na(cell_type)
  n_samples_celltypes_available <- length(unique(model_input_df$SUBJID))
  model_input_df$n_samples_celltypes_available <- n_samples_celltypes_available
  model_input_df$frac_celltypes_available <- n_samples_celltypes_available/n_samples_total
  
  # subset to pos_feature for which we found significant donor age association
  # t_pos <- data_to_test_df[which(data_to_test_df$tissue == tissue & data_to_test_df$BB_E_bonferroni <= 0.05), "Pos"]
  t_pos <- data_to_test_df[which(data_to_test_df$tissue == tissue), "Pos"]
  t_pos <- t_pos %>% pull(Pos)
  #model_input_df <- model_input_df[which(model_input_df$Pos %in% t_pos), ]
  
  # apply test below to each eQTL position feature pair
  p_res_list <- list()
  for(p in seq_along(t_pos)){
    # subset for pos_feature pairs
    position <- t_pos[p]
    p_df <- model_input_df[which(model_input_df$Pos %in% position), ] 
    
    # apply test per cell_type
    cell_types <- unique(p_df$cell_type) 
    c_res_list <- list()
    for(c in seq_along(cell_types)){ 
      # subset per cell_type
      cell_type <- cell_types[c]
      c_df <- p_df[which(p_df$cell_type == cell_type), ]
      # init empty result if errors would occur
      na_df <- data.frame(tissue = tissue,
                          Pos = position,
                          cell_type = cell_type,                 
                          n_samples = length(unique(c_df$SUBJID)),
                          n_donors_cell_type_file = length(unique(c_df$SUBJID)),
                          mean_cell_type_proportion = NA,
                          frac_celltypes_available = unique(c_df$frac_celltypes_available),
                          pass_testing_filters = FALSE,
                          lm_m1_interaction_estimate = NA,
                          lm_m1_interaction_pval = NA,
                          LR_pval = NA)  
      # only run test if variance is > 0 (i.e. if there is all 0 % for one c)
      # let's require 5 donors per cell type having more than 0 to be considered 
      if(var(c_df[ ,"proportion"]) > 0.00 && nrow(c_df[which(c_df$proportion == 0), ]) <= 5){
      ### define model formulas
        single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
        if(unique(p_df$tissue) %in% single_sex_tissues){   
          lm_m0_formula <- paste0("sum_heteroplasmic_level ~ AGE + proportion + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
          lm_m1_formula <- paste0("sum_heteroplasmic_level ~ AGE + proportion + AGE : proportion + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        } else {
          lm_m0_formula <- paste0("sum_heteroplasmic_level ~ AGE + proportion + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
          lm_m1_formula <- paste0("sum_heteroplasmic_level ~ AGE + proportion + AGE : proportion + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
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
        # Use tryCatch to attempt to create the desired data frame
        c_res_df <- tryCatch({
          data.frame(tissue = tissue,
                     Pos = position,
                     cell_type = cell_type,                 
                     n_samples = length(unique(c_df$SUBJID)),
                     n_donors_cell_type_file = length(unique(c_df$SUBJID)),
                     mean_cell_type_proportion = mean(c_df$proportion),
                     frac_celltypes_available = unique(c_df$frac_celltypes_available),
                     pass_testing_filters = TRUE,
                     lm_m1_interaction_estimate = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("AGE:proportion")), "Estimate"],
                     lm_m1_interaction_std_error = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("AGE:proportion")), "Std. Error"],
                     lm_m1_interaction_pval = res_m1$coefficients[which(rownames(res_m1$coefficients) == paste0("AGE:proportion")), "Pr(>|t|)"],
                     LR_pval = LR_test[2,'Pr(>Chisq)'])
        }, error = function(e) {
          # If an error occurs, return the NA-filled data frame
          print("An error occurred, returning NA-filled data frame.")
          print(e)
          return(na_df)
        })
        
      } else {
        # return NAs for those cell types not present or tested
        c_res_df <- na_df
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

#+ do study-wide bonferroni correction ----------------------------------------
res_df <- res_df %>%
  filter(pass_testing_filters == TRUE) %>%
  mutate(LM_A_bonferroni =  p.adjust(LR_pval, method = 'bonferroni')) 

#+ print significant summaries to console --------------------------------------
# significant across all tests
res_df %>% 
  filter(LM_A_bonferroni <= 0.05) %>%
  kbl() %>%
  kable_paper(bootstrap_options = "striped", full_width = F)
# signifcant in previous donor age findings
left_join(res_df,donor_age_sig_df %>% select(tissue,Pos) %>% mutate(donor_age_sig = TRUE)) %>% filter(LM_A_bonferroni <= 0.05)


#+ save to disc ----------------------------------------------------------------
write_tsv(res_df,paste0(save_path,"2024_08_14_LM_A_donor_age_all_tissues_cell_type_deconvolution.tsv"))
