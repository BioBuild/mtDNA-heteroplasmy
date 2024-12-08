# author:   simon.wengert@gmail.com
source("~/scripts/utils/global_settings.R")
library("ggpubr")

#+ load data -------------------------------------------------------------------
# the cleaned set of heteroplasmies 
gt_rna_var_df <- readRDS(paste0(g_ops$metadata_dir,"2024_04_28_gt_rna_var_annotated.rds"))
n_samples_df <- gt_rna_var_df %>% 
  distinct(tissue,biospecimen_repository_sample_id) %>% 
  count(tissue) %>% 
  mutate(n_samples = n)
gt_rna_var_df <- left_join(gt_rna_var_df,n_samples_df)
# the genotype files for all tissues, which also contains all the heteroplasmy
# variation detected in the data. 
tissues <- unique(gt_rna_var_df$tissue)
all_t_list <- list()
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("load and concatenate heteroplasmy genotype files per tissue #",t,"/",length(tissues)))
  t_genotypes_df <- read_tsv(paste0("~/data/mtDNA_variants/metadata/genotype_files/split_by_tissue/",tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)
  t_genotypes_df$tissue <- tissue
  all_t_list[[t]] <- t_genotypes_df
}
gt_genotypes_df <- bind_rows(all_t_list)
# gt_genotypes_df %>% distinct(tissue,Pos) %>% nrow(.)
# [1] 4334 --> number mentioned in the main text
# quick sanity check if the number of positions matches our expectation for the
# data we expect to load.
# length(unique(gt_genotypes_df$Pos))
# [1] 556 --> it does
# get RNA modifications for plotting
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))

################################################################################
###       TEST FOR HETEROSKEDASTICITY HET LEVEL (VAF) VS COVERAGE            ###              
################################################################################
#+ apply breusch pagan test see if there is any evidence for oversdipersion ----
filename_overdisp_table <-"2024_10_24_gtex_BP_test_for_overdispersion_all_pos_as_used_in_association_testing.csv"
if(!file.exists(filename_overdisp_table)){
  # get data as used in association tests 
  all_t_df <- gt_genotypes_df %>%    
      dplyr::select(tissue,SUBJID,AGE,Pos,Coverage,sum_heteroplasmic_level) %>%
      distinct()
  res_list <- list()
  tissues <- unique(all_t_df$tissue)
  for(t in seq_along(tissues)){
    cat(paste0("tissue is #",t,"/",length(tissues)," --------------- \n"))
    tissue <- tissues[t]
    t_df <- all_t_df[which(all_t_df$tissue == tissue), ] 
    t_res_list <- list()
    positions <- unique(t_df$Pos)
    for(p in seq_along(positions)){
      # subset for position
      position <- positions[p]
      p_df <- t_df[which(t_df$Pos == position), ] 
      # fit LM per position
      model_formula <- "sum_heteroplasmic_level ~ Coverage"
      # lm_m1 <- lm(sum_heteroplasmic_level ~ AGE, data = p_df)
      # 2024-10-24: we talk about heteroplasmy ~ Coverage throughout in the 
      # manuscript!!!
      lm_m1 <- lm(as.formula(model_formula), data = p_df)
      # linear model testing for overdispersion
      lm_bptest <- lmtest::bptest(lm_m1)
      lm_bptest_method <- lm_bptest$method
      lm_bptest_statistic <- lm_bptest$statistic
      lm_bptest_p_val <- lm_bptest$p.value
      # assemble result
      pos_res_df <- tibble(tissue = tissue, 
                           Pos = position,
                           LM_formula = model_formula,
                           method = lm_bptest_method,
                           LM_A_bptest_p_val = lm_bptest_p_val,
                           LM_A_bptest_statistic = lm_bptest_statistic)
      # return result
      t_res_list[[p]] <- pos_res_df
    }
    t_res_df <- bind_rows(t_res_list)  
    # tissue wide multiple testing correction
    t_res_df$LM_A_bptest_bonferroni_tissue_wide <- p.adjust(t_res_df$LM_A_bptest_p_val, method = "bonferroni")
    res_list[[t]] <- t_res_df
  }
  res_df <- bind_rows(res_list)
  # save to disc
  write_csv(res_df,filename_overdisp_table)
}else{
  res_df <- read_csv(filename_overdisp_table)
}
# multiple testing correction
res_df$LM_A_bptest_bonferroni_study_wide <- p.adjust(res_df$LM_A_bptest_p_val, method = "bonferroni")
# CHOOSE WHICH BONFERONNI TO GO WITH
res_df$LM_A_bptest_bonferroni <- res_df$LM_A_bptest_bonferroni_tissue_wide
# go on with the script + files for supplement
n_sig <- res_df %>% filter(LM_A_bptest_bonferroni <= 0.05) %>% nrow(.)
n_tests <- nrow(res_df)
# make a histogram of BP test p-values for all tissues
res_df$plot_tissue <- paste0("all tissues (# tests: ",n_tests,", # bonferroni <= 0.05: ",n_sig,")")
# section out significant ones, rename to human readable columns and cleanup
res_sig_df <- res_df %>% 
  filter(LM_A_bptest_bonferroni <= 0.05) %>%
  dplyr::select(tissue,Pos,method,
                LM_formula,
                test_statistic = LM_A_bptest_statistic,
                p_value = LM_A_bptest_p_val,
                bonferroni = LM_A_bptest_bonferroni)
# save signigicant BP test results to disc for supp table 
write_tsv(res_sig_df,"2024_10_24_gtex_BP_test_for_overdispersion_bonferroni_passing.tsv")

################################################################################
###       EXTRACT REAL DATA EQUIVLALENTS OF SIM PARAMS ACROSS ALL SITES      ###              
################################################################################

data_type <- "genotypes"
# data_type <- "heteroplasmy_pos"
if(data_type == "heteroplasmy_pos"){
  assess_df <- gt_rna_var_df
} else if (data_type == "genotypes"){
  assess_df <- gt_genotypes_df %>% 
    drop_na()
}
# loop thorough tissue-pos to find linear trends if there are any:
tissues <- unique(assess_df$tissue)  
res_list <- list()
all_pos_res_list <- list()
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("tissue is #",t,"/",length(tissues)," --------------- \n"))
  tissue_df <- assess_df[which(assess_df$tissue == tissue), ]
  positions <- unique(tissue_df$Pos)
  # go to loop 
  t_res_list <- list()
  for(p in seq_along(positions)){
    position <- positions[p]
    p_df <- tissue_df[which(tissue_df$Pos == position), ] 
    
    # we are going to test for relationship between 
    ## (i) coverage and heteroplasmy and ...
    cov_het_linear_model <- lm(Coverage ~ sum_heteroplasmic_level, data = p_df)
    # get summary 
    cov_het_lm <- summary(cov_het_linear_model)
    # get estimate
    cov_het_pos_estimate <- tryCatch({ 
      estimate <- cov_het_lm$coefficients[2, 1] 
      ifelse(is.na(estimate), NA, estimate)
    }, error = function(e) {NA})
    # get standard error
    cov_het_pos_se <- tryCatch({ 
      se <- cov_het_lm$coefficients[2, 2]
      ifelse(is.na(se), NA, se)
    }, error = function(e) {NA})
    # get p-value
    cov_het_pos_pval <- tryCatch({ 
      p_value <- cov_het_lm$coefficients[2, 4] 
      ifelse(is.na(p_value), NA, p_value)
    }, error = function(e) {NA})
    
    ## (ii) coverage vs. age
    # get slope coefficient of age vs coverage
    cov_age_linear_model <- lm(Coverage ~ AGE, data = p_df)
    # get summary 
    cov_age_lm <- summary(cov_age_linear_model)
    # get estimate
    cov_age_pos_estimate <- tryCatch({ 
      estimate <- cov_age_lm$coefficients[2, 1] 
      ifelse(is.na(estimate), NA, estimate)
      }, error = function(e) {NA})
    # get standard error
    cov_age_pos_se <- tryCatch({ 
      se <- cov_age_lm$coefficients[2, 2]
      ifelse(is.na(se), NA, se)
      }, error = function(e) {NA})
    # get p-value
    cov_age_pos_pval <- tryCatch({ 
      p_value <- cov_age_lm$coefficients[2, 4] 
      ifelse(is.na(p_value), NA, p_value)
      }, error = function(e) {NA})
    # get mean coverage at this position across donors
    mean_cov <- mean(as.integer(p_df$Coverage))
    # get mean heteroplasmy
    mean_het <- mean(p_df$sum_heteroplasmic_level)
    # get var heteroplasmy
    var_het <- var(p_df$sum_heteroplasmic_level)
    # build result
    tmp_df <- tibble(Pos = position,
                     tissue = unique(p_df$tissue),
                     lm_cov_het_estimate = cov_het_pos_estimate,
                     lm_cov_het_std_error = cov_het_pos_se,
                     lm_cov_het_p_value = cov_het_pos_pval, 
                     lm_cov_age_estimate = cov_age_pos_estimate,
                     lm_cov_age_std_error = cov_age_pos_se,
                     lm_cov_age_p_value = cov_age_pos_pval,
                     mean_coverage = mean_cov,
                     r_squared = cov_age_lm$r.squared,
                     intercept = coef(cov_age_linear_model)[1], # base_coverage is just the intercept
                     slope = coef(cov_age_linear_model)[2],
                     mean_het = mean_het,
                     var_het = var_het)
    # return result
    t_res_list[[p]] <- tmp_df
  }
  # collapse list to data frame
  t_res_df <- bind_rows(t_res_list) 
  # let's do multiple testing correction here
  ## coverage heteroplasmy relationship
  t_res_df$lm_cov_het_fdr <- p.adjust(t_res_df$lm_cov_het_p_value, method = "fdr")
  t_res_df <- t_res_df %>% 
    mutate(cov_het_sig = case_when(lm_cov_het_fdr <= 0.05 ~ TRUE,TRUE ~ FALSE))
    ## coverage age relationship
  t_res_df$lm_cov_age_fdr <- p.adjust(t_res_df$lm_cov_age_p_value, method = "fdr")
  t_res_df <- t_res_df %>% 
    mutate(cov_age_sig = case_when(lm_cov_age_fdr <= 0.05 ~ TRUE,TRUE ~ FALSE))
  # only keep significant linear relationships between age ~ coverage
  # since this is the model we assume for simulations
  t_cov_age_res_df <- t_res_df %>% filter(cov_age_sig == TRUE)
  # get ranges real data params for doing simulations later
  t_slope_summary <- summary(t_cov_age_res_df$slope)
  t_mean_het_summary <- summary(t_cov_age_res_df$mean_het)
  t_var_het_summary <- summary(t_cov_age_res_df$var_het)
  # build res df and return
  t_tmp_df <- data.frame(tissue = tissue,
                         slope_min = as.double(t_slope_summary[1]),
                         slope_median = as.double(t_slope_summary[3]),
                         slope_max = as.double(t_slope_summary[6]),
                         mean_het_min = as.double(t_mean_het_summary[1]),
                         mean_het_median = as.double(t_mean_het_summary[3]),
                         mean_het_max = as.double(t_mean_het_summary[6]),
                         var_het_min = as.double(t_var_het_summary[1]),
                         var_het_median = as.double(t_var_het_summary[3]),
                         var_het_max = as.double(t_var_het_summary[6]))
  res_list[[t]] <- t_tmp_df
  # keep all pos for debugging underlying positions
  all_pos_res_list[[t]] <- t_res_df
}
# collaps to data frame
res_df <- bind_rows(res_list)
all_pos_res_df <- bind_rows(all_pos_res_list)
# save this table
write_csv(res_df,paste0("2024_09_26_heteroplasmy_param_ranges_real_tissues_basis_for_simulation_from_",data_type,".csv"))
write_csv(all_pos_res_df,paste0("2024_09_26_gtex_all_pos_slope_and_r_sqr_observed_from_",data_type,".csv"))
# **NOTE:** if there is only NAs for a given tissue this means there are no FDR
# significant age-coverage trends found for this tissue. (in res_df in particular)

# subset supp table for COV ~ HET
cov_het_supp_table_df <- all_pos_res_df %>%
  filter(cov_het_sig == TRUE) %>%
  select(Pos,tissue,
         estimate = lm_cov_het_estimate,
         std_error = lm_cov_het_std_error,
         p_value = lm_cov_het_p_value,
         bonferroni = lm_cov_het_fdr)
# save to disc 
write_tsv(cov_het_supp_table_df,"2024_07_05_gtex_cov_het_supp_table_bonferroni_passing.tsv")
# 
# subset supp table for COV ~ AGE
cov_age_supp_table_df <- all_pos_res_df %>%
  filter(cov_age_sig == TRUE) %>%
  select(Pos,tissue,
         estimate = lm_cov_age_estimate,
         std_error = lm_cov_age_std_error,
         p_value = lm_cov_age_p_value,
         bonferroni = lm_cov_age_fdr)
# save to disc 
write_tsv(cov_age_supp_table_df,"2024_07_05_gtex_cov_age_supp_table_bonferroni_passing.tsv")

################################################################################
###            VISUALISE REAL PARAMS USING MIXER STYLE BOXPLOTS              ###              
################################################################################

sim_params_df <- read_tsv("2024_05_01_sim_param_combinations_used_for_model_evaluation.tsv",show_col_types = FALSE)
# pull out extreme values we set for sim params for highlighting below
slope_min <- min(sim_params_df$slope)
slope_max <- max(sim_params_df$slope)
r_sqr_min <- min(sim_params_df$r_sqr)
# make same number of digits in the plot
#r_sqr_max <- sprintf(max(sim_params_df$r_sqr), fmt = '%#.2f')
r_sqr_max <- max(sim_params_df$r_sqr)
mean_coverage_min <- min(sim_params_df$mean_coverage)
mean_coverage_max <- max(sim_params_df$mean_coverage)
mean_beta_min <- min(sim_params_df$mean_beta)
# make same number of digits in the plot
#mean_beta_max <- sprintf(max(sim_params_df$mean_beta), fmt = '%#.2f')
mean_beta_max <- max(sim_params_df$mean_beta)


# now let's make a boxplot for each parameter we care about, summary used for simulations
# a) subplot for value ranges observed for slope:
min_max_data <- filtered_data %>% filter(slope == max(slope) | slope == min(slope))
# let's keep a record of extreme values for reference/supplement
min_max_data$extreme_value_type <- "slope"
min_max_data$sim_params <- c(slope_min,slope_max)
# make plot
p_1 <- filtered_data %>%
  mutate(is_extreme = case_when(Pos %in% min_max_data$Pos & tissue %in% min_max_data$tissue ~ TRUE,
                                TRUE ~ FALSE)) %>%
  ggplot(aes(x = "", y = slope)) +
  geom_hline(yintercept = slope_min, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_hline(yintercept = slope_max, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_boxplot(outlier.shape = NA) +
  # use 2 geom point layers for making sure red dots are on top of black ones.
  geom_point(fill = "black", colour = "white", stroke = 0.2, pch = 21, size = 2.5, 
             position = position_jitter(width = 0.15)) +  # Black points
  geom_point(data = . %>% filter(is_extreme == TRUE),fill = "red", colour = "white", stroke = 0.2, pch = 21, size = 3.5, 
  ) + 
  geom_text(data = min_max_data,
            aes(label = sim_params, y = sim_params),
            vjust = -2.6,
            hjust = 1.15,
            color = "darkblue", fontface = "bold", size = 10/.pt) +
  ylim(-420,390) +
  labs(y = "", x = "slope") +
  coord_flip() +
  p_ops$grid_theme +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0.5, 2, 0, 2, "pt"),
        plot.= element_blank()) 
# b) subplot for value ranges observed for r_sqr:
min_max_data <- filtered_data %>% filter(r_squared == max(r_squared) | r_squared == min(r_squared))
# let's keep a record of extreme values for reference/supplement
min_max_data$extreme_value_type <- "r_sqr"
min_max_data$sim_params <- c(r_sqr_min,r_sqr_max)
# make plot
p_2 <- filtered_data %>%
  mutate(is_extreme = case_when(Pos %in% min_max_data$Pos & tissue %in% min_max_data$tissue ~ TRUE,
                                TRUE ~ FALSE)) %>%
  ggplot(aes(x = "", y = r_squared)) +
  geom_hline(yintercept = r_sqr_min, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_hline(yintercept = r_sqr_max, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_boxplot(outlier.shape = NA) +
  # use 2 geom point layers for making sure red dots are on top of black ones.
  geom_point(fill = "black", colour = "white", stroke = 0.2, pch = 21, size = 2.5, 
             position = position_jitter(width = 0.15)) +  # Black points
  geom_point(data = . %>% filter(is_extreme == TRUE),fill = "red", colour = "white", stroke = 0.2, pch = 21, size = 3.5, 
  ) + 
  geom_text(data = min_max_data,
            # below is for consistently plotting 2 digits
            aes(label = sprintf("%0.2f",sim_params), y = sim_params),
            vjust = -2.6,
            hjust = 1.15,
            color = "darkblue", fontface = "bold", size = 10/.pt) +
  ylim(0.005,0.1175) +
  labs(y = "",x = "r_sqr") +
  coord_flip() +
  p_ops$grid_theme +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0.3, 2, 0, 2, "pt")) 
# c) subplot for value ranges observed for mean_coverage:
min_max_data <- filtered_heteroplasmies_df %>% filter(mean_coverage == max(mean_coverage) | mean_coverage == min(mean_coverage))
# let's keep a record of extreme values for reference/supplement
min_max_data <- min_max_data %>% dplyr::distinct(Pos,tissue,mean_coverage)
min_max_data$extreme_value_type <- "mean_coverage"
min_max_data$sim_params <- c(mean_coverage_min,mean_coverage_max)
# make plot
p_3 <- filtered_heteroplasmies_df %>%
  distinct(tissue,Pos,mean_coverage) %>%
  mutate(is_extreme = case_when(Pos %in% min_max_data$Pos & tissue %in% min_max_data$tissue ~ TRUE,
                                TRUE ~ FALSE)) %>%
  ggplot(aes(x = "", y = mean_coverage)) +
  geom_hline(yintercept = mean_coverage_min, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_hline(yintercept = mean_coverage_max, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_boxplot(outlier.shape = NA) +
  # use 2 geom point layers for making sure red dots are on top of black ones.
  geom_point(fill = "black", colour = "white", stroke = 0.2, pch = 21, size = 2.5, 
             position = position_jitter(width = 0.15)) +  # Black points
  geom_point(data = . %>% filter(is_extreme == TRUE),fill = "red", colour = "white", stroke = 0.2, pch = 21, size = 3.5, 
  ) + 
  geom_text(data = min_max_data,
            aes(label = sim_params, y = sim_params),
            vjust = -2.6,
            hjust = 1.2,
            color = "darkblue", fontface = "bold", size = 10/.pt) +
  labs(y = " ", x = "mean_coverage") +
  ylim(-200,20000) +
  coord_flip() +
  p_ops$grid_theme +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0.3, 2, 0, 2, "pt")) 
# d) subplot for value ranges observed for mean_heteroplasmy:
min_max_data <- filtered_heteroplasmies_df %>% filter(mean_het == max(mean_het) | mean_het == min(mean_het))
# let's keep a record of extreme values for reference/supplement
min_max_data <- min_max_data %>% dplyr::distinct(Pos,tissue,mean_het)
min_max_data$extreme_value_type <- "mean_het"
min_max_data$sim_params <- c(mean_beta_min,mean_beta_max)
# make plot
p_4 <- filtered_heteroplasmies_df %>%
  distinct(tissue,Pos,mean_het) %>%
  mutate(is_extreme = case_when(Pos %in% min_max_data$Pos & tissue %in% min_max_data$tissue ~ TRUE,
                                TRUE ~ FALSE)) %>%
  ggplot(aes(x = "", y = mean_het)) +
  geom_hline(yintercept = mean_beta_min, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_hline(yintercept = mean_beta_max, linetype = "dashed", linewidth = 0.8, colour = "darkblue") +
  geom_boxplot(outlier.shape = NA) +
  # use 2 geom point layers for making sure red dots are on top of black ones.
  geom_point(fill = "black", colour = "white", stroke = 0.2, pch = 21, size = 2.5, 
             position = position_jitter(width = 0.15)) +  # Black points
  geom_point(data = . %>% filter(is_extreme == TRUE),fill = "red", colour = "white", stroke = 0.2, pch = 21, size = 3.5, 
  ) + 
  geom_text(data = min_max_data,
            # below is for consistently plotting 2 digits
            aes(label = sprintf("%0.2f",sim_params), y = sim_params),
            vjust = -2.6,
            hjust = 1.12,
            color = "darkblue", fontface = "bold", size = 10/.pt) +
  ylim(-0.001,0.505) +
  labs(y = " ", x = "mean_het") +  
  coord_flip() +
  p_ops$grid_theme +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(0.3, 2, 0, 2, "pt")) 
# e) assemble the plots into panel
p_panel_data_observed_ranges_for_simulations <- cowplot::plot_grid(p_1,p_2,p_3,p_4,ncol = 1)
# save plots to disc
ggsave(p_panel_data_observed_ranges_for_simulations, 
       filename = paste0("~/git/mtDNA_variants/paper/figure_plots/simulations/2024_07_03_p_panel_data_observed_ranges_for_simulations_with_sim_param_boundaries.png"),
       width = 6,
       height = 5.5,
       dpi = "retina")   
ggsave(p_panel_data_observed_ranges_for_simulations, 
       filename = paste0("~/git/mtDNA_variants/paper/figure_plots/simulations/2024_07_03_p_panel_data_observed_ranges_for_simulations_with_sim_param_boundaries.pdf"),
       width = 5.8,
       height = 5.5,
       useDingbats = FALSE,
       dpi = 300)   
