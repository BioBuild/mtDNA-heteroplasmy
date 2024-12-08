# author:  simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library(lmtest)

model_input_path <- "/split_by_tissue/"
model_output_path <- "/donor_age/"
tissues <- list.dirs(model_input_path, full.names = FALSE)[-1]

#+ loop through tissues and fit lm model per tissue to find candidate SNPs -----
tissue_res_list <- list()
for(t in seq_along(tissues)){
tissue <- tissues[t]
cat(paste0("\n#tissues is ",t," out of #",length(tissues),"\n"))
model_input_data <- read_tsv(paste0(model_input_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)

#+ drop NAs for now doing complete case analysis -------------------------------
model_input_data <- model_input_data %>% drop_na()


#+ filter out positions 0 variance ---------------------------------------------
model_input_data  <- 
  model_input_data %>%
    filter(var_per_pos > f_ops$site_variance_larger_than)


#+ filter for minimum number of donors -----------------------------------------
min_donors_pos <- 
  model_input_data  %>%
      group_by(Pos) %>%
      count(Pos) %>% 
        filter(n >= f_ops$site_n_donors_min) %>%
      ungroup() %>%  
      pull(Pos)

    # aplly LM to test for donor age relationship 
    # heteroplasmy candidates
    positions <- unique(model_input_data$Pos)
    positions <- positions[which(positions %in% min_donors_pos)]
    # fitting the model
    res_list <- list()
    for (p in seq_along(positions)){
        position <- positions[p]
        pos_input_data <- 
            model_input_data %>%
                filter(Pos == position)
        lm_m0 <- lm(sum_heteroplasmic_level ~ SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_data)
        lm_m1 <- lm(sum_heteroplasmic_level ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_data)
        lm_m1_summary_df <- summary(lm_m1)
        LR_test <- lrtest(lm_m1,lm_m0) 
        tmp_df <- data.frame(tissue = tissue,
                             Pos = position,
                             estimate = lm_m1$coefficients["AGE"],
                             std_error = lm_m1_summary_df$coefficients["AGE","Std. Error"],
                             analytical_p_value =  LR_test[2,  'Pr(>Chisq)'],
                             n_donors = length(unique(pos_input_data$SUBJID)),
                             n_het = nrow(pos_input_data[which(pos_input_data$signal_type == "heteroplasmy"), ]),
                             n_non_het = nrow(pos_input_data[which(pos_input_data$signal_type == "non_heteroplasmy"), ]),
                             n_missing = nrow(pos_input_data[which(pos_input_data$signal_type == "missing_value"), ]))
        res_list[[p]] <- tmp_df
    }
    # collapse to data frame
    lm_t_df <- bind_rows(res_list)
    # apply multiple testing correction -- see below, we now do this study-wide
    #lm_t_df$analytical_FDR <-  p.adjust(lm_t_df$analytical_p_value, method = 'fdr')
    #lm_t_df$analytical_bonferroni <- p.adjust(lm_t_df$analytical_p_value, method = 'bonferroni')
    # store in list
    tissue_res_list[[t]] <- lm_t_df
}
tissue_res_df <- bind_rows(tissue_res_list)
tissue_res_df$analytical_FDR <-  p.adjust(tissue_res_df$analytical_p_value, method = 'fdr')
tissue_res_df$analytical_bonferroni <- p.adjust(tissue_res_df$analytical_p_value, method = 'bonferroni')


#+ save results to disc --------------------------------------------------------
write_tsv(tissue_res_df,paste0(model_output_path,"2024_07_23_gt_LM_A_donor_age_all_tissues_common_het_pos.tsv"))
