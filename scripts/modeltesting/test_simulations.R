# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library("tidyverse")
library("ggpubr") # for table plots of simulation parameters
library("ggsignif") # for significance bars

# set notation to full numnbers for plotting
options(scipen = 999)

#+ load data -------------------------------------------------------------------
all_sim_df <- read_tsv("2024_05_01_merged_simulated_data.tsv")
all_sim_df <-
  all_sim_df %>%
    # make coverage full integers
    mutate(coverage = round(coverage)) %>%
    mutate(AD = het_counts,
           BP = coverage - het_counts)
# make data frame for subsetting using base R below
all_sim_df <- as.data.frame(all_sim_df)
# this is simulation parameters for all the simulations succ. run
sim_params_df <- 
    all_sim_df %>%
        distinct(row_id,date,change_param,step_param,n_donors,mean_age,sd_age,min_age_bound,
                 max_age_bound,min_coverage,r_sqr,slope,mean_coverage,mean_beta,
                 var_beta)
# save these parameters of models we actually fitted to table for supplement
write_tsv(sim_params_df,"2024_05_01_sim_param_combinations_used_for_model_evaluation.tsv")
if(F){
# now adding the in which parameter was intended to be changed (used for showing 
# a suitable range for each and also finding which range makes sense).
#sim_param_specs_df <- read_tsv("~/git/mtDNA_variants/metadata/utils/het_sim_params_df_2.txt")
sim_param_specs_df <- 
  read_csv("~/git/mtDNA_variants/metadata/utils/het_sim_params.csv") %>%
  filter(run == TRUE)
sim_params_df <- 
  left_join(sim_params_df,
            sim_param_specs_df %>%
              select(sim_param_row_id = row_id,
                     sim_param_base_coverage = base_coverage,
                     change_parameter = change_param),
            by = c("sim_param_row_id","sim_param_base_coverage"))
}


################################################################################
###         fit LM & BB models to simulated data  - analytical pvals         ###
################################################################################
if(!file.exists(all_models_path)){
  # get iterator for params of interest
  iterate_rows <- sim_params_df %>% pull(row_id)
  all_models_list <- list()
  for(i in seq_along(iterate_rows)){
    # subset simulated data to one single parameter combination
    cat(paste0("\n fitting models to simulated data fro row_id: #",i,"/",length(iterate_rows)))
    sim_df <- all_sim_df %>% filter(row_id == iterate_rows[i])
    # check if there's heteroplasmy for this set of params at all
    if(length(unique(sim_df$het_level)) == 1 && unique(sim_df$het_level) == 0){
      next
    }else{
      #+ get analytical p-values for assessing model calibration ---------------
      # fit LM and BB per position
      positions <- unique(sim_df$Pos)
      res_list <- list()
      for(p in seq_along(positions)){
        # subset for position
        position <- positions[p]
        p_df <- sim_df %>% filter(Pos == position)
        # fit LM per position
        lm_m0 <- lm(het_level ~ 1, data = p_df)
        # LM: fit test_reltionship 
        if(test_relationship == "het_level ~ coverage"){
           lm_m1 <- lm(het_level ~ coverage, data = p_df)
        } else if(test_relationship == "het_level ~ age") {
           lm_m1 <- lm(het_level ~ age, data = p_df)
        } else if(test_relationship == "het_level ~ age + coverage") {
           lm_m1 <- lm(het_level ~ age + coverage, data = p_df)
        } else {
          stop("LM: please provide a valid model formula")
        }
        lm_estimate <- as.double(lm_m1$coefficients[2])
        # get LRs: lmtest is still not working. So let's derive LR manually
        logLL_lm_m0 <- logLik(lm_m0)
        logLL_lm_m1 <- logLik(lm_m1)
        # Calculate the likelihood ratio test statistic
        LR_test_statistic <- 2*(logLL_lm_m1 - logLL_lm_m0)
        # Calculate degrees of freedom
        df <- attr(logLL_lm_m1, "df") - attr(logLL_lm_m0, "df")
        # Perform the likelihood ratio test
        lm_p_val <- as.double(1 - pchisq(LR_test_statistic, df))
        # linear model testing for overdispersion
        lm_bptest <- lmtest::bptest(lm_m1)
        lm_bptest_p_val <- lm_bptest$p.value
        # fit BB_per position
        bb_m0 <- aod::betabin(cbind(AD,BP) ~ 1, ~ 1, data = p_df)
        # BB: fit test_reltionship 
        if(test_relationship == "het_level ~ coverage"){
           bb_m1 <- aod::betabin(cbind(AD,BP) ~ coverage, ~ 1, data = p_df)
        } else if(test_relationship == "het_level ~ age") {
           bb_m1 <- aod::betabin(cbind(AD,BP) ~ age, ~ 1, data = p_df)
        } else if(test_relationship == "het_level ~ age + coverage") {
           bb_m1 <- aod::betabin(cbind(AD,BP) ~ age + coverage, ~ 1, data = p_df)
        } else {
          stop("BB: please provide a valid model formula")
        }
        # get effect size
        bb_tmp <- show(bb_m1)
        bb_estimate <- bb_tmp@Coef[2, "Estimate"]
        # get LL ratios
        LR_bb = 2*(bb_m1@logL - bb_m0@logL)
        # get pvalue
        bb_p_val <-  pchisq(LR_bb, df = 1, lower.tail = FALSE, log.p = FALSE) 
        # build res_df
        pos_res_df <- tibble(Pos = position,
                             test_relationship = test_relationship,
                             value_type = "p_analytical",
                             LM_estimate = lm_estimate,
                             BB_estimate = bb_estimate,
                             LM_A_p_val = lm_p_val,
                             LM_A_bptest_p_val = lm_bptest_p_val,
                             BB_A_p_val = bb_p_val)
        # return result
        res_list[[p]] <- pos_res_df
      }
    }
    res_df <- bind_rows(res_list)
    # let's do a study-wide bonferroni
    res_df$LM_A_bonferroni <- p.adjust(res_df$LM_A_p_val, method = "bonferroni")
    res_df$BB_A_bonferroni <- p.adjust(res_df$BB_A_p_val, method = "bonferroni")
    # keep track of false positive rates
    res_df$LM_A_p_vals_small <- res_df %>% filter(LM_A_p_val <= 0.05) %>% nrow(.)
    res_df$BB_A_p_vals_small <- res_df %>% filter(BB_A_p_val <= 0.05) %>% nrow(.)
    # keep track of simulation paramters for join ops below
    res_df$row_id <- unique(sim_df$row_id)
    # stash results for later
    all_models_list[[i]] <- res_df
  }
  all_models_df <- bind_rows(all_models_list)
  # save all models df
  write_tsv(all_models_df,all_models_path)
  # tell me if run, saved and done  
  cat(paste0("\n fitting models to simulated data SAVED AND DONE ________________"))
} else {
  all_models_df <- read_tsv(all_models_path)
}

# also annotate model outputs with change parameters
all_models_df <-
  left_join(all_models_df,
            sim_params_df,
            by = c("row_id"))

#+ make correltion plot of FP rates for supplement -----------------------------
# make plot
p_scatter_plot_fp_rates_lm_and_bb_models <- 
  plot_delta_fp_rate_df %>%
    ggplot(aes(BB_fp_rate,LM_fp_rate,label = simulation_run)) +
    geom_point(shape = 21,stroke = 2,size = 2.5) +
    geom_label_repel(nudge_x = 0.002,nudge_y = 0.001) +  
    geom_abline(lty = "dashed") +
    annotate("segment", x = 0.025, xend = 0.05, y = 0.05, yend = 0.05, 
             colour = "darkred", linewidth = 1.5, linetype = "dashed") +
    annotate("segment", y = 0.025, yend = 0.05, x = 0.05, xend = 0.05, 
             colour = "darkred", linewidth = 1.5, linetype = "dashed") +
    xlim(0.025,0.075) +
    ylim(0.025,0.075) +
    labs(title = "Comparing FP rates LM and BB",
         x = "BB FP rate",
         y = "LM FP rate") +
    p_ops$my_theme
# save to disc 
ggsave(p_scatter_plot_fp_rates_lm_and_bb_models, 
       filename = paste0("2024_07_02_p_scatter_plot_fp_rates_lm_and_bb_models.png"),
       width = 6,
       height = 5.5,
       dpi = 300) 
ggsave(p_scatter_plot_fp_rates_lm_and_bb_models, 
       filename = paste0("2024_07_02_p_scatter_plot_fp_rates_lm_and_bb_models.pdf"),
       width = 6,
       height = 5.5,
       useDingbats = FALSE,
       dpi = 300) 
