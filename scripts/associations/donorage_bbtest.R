#' author: 'Simon Wengert (simon.wengert@helmholtz-munich.de)'
source("~/scripts/utils/global_settings.R")
source('~/scripts/associations/function_beta_binomial_empirical_p_values.R')
library("data.table")

metadata_dir <- "/split_by_tissue/"
lm_res_all_tissues_df <- read_tsv("2024_07_23_gt_LM_A_donor_age_all_tissues_common_het_pos.tsv", show_col_types = FALSE)

#+ loop through tissues --------------------------------------------------------
tissues <- list.dirs(g_ops$metadata_dir, full.names = FALSE)[-1]
t_res_list <- list()
for(t in seq_along(tissues)){
    # set tissue
    tissue <- tissues[t]
    cat(paste0("\nStarting BB_A donor age for tissue #",t,"/",length(tissues),": ",tissue))
    # let's keep an eye on runtime
    start_time <- Sys.time()
    # load tissue genotypes 
    t_genotypes <- read_tsv(paste0(g_ops$metadata_dir,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)
    #
    # subset to positions with LM analtical significance
    # pos_idx <- which(lm_res_all_tissues_df$tissue == tissue & lm_res_all_tissues_df$analytical_FDR <= 0.05)
    pos_idx <- which(lm_res_all_tissues_df$tissue == tissue & lm_res_all_tissues_df$analytical_bonferroni <= 0.05)
    positions <- as.data.frame(lm_res_all_tissues_df)[pos_idx, "Pos"]
    #
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
    # only fit model for positions > 0
    if(length(positions) > 0){
      # define association test formula depending if tissue is present in 
      # one or both of the sexes
      single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
      if(unique(input_data$tissue) %in% single_sex_tissues){   
        null_model_formula <- c("cbind(AD, BP) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        alt_model_formula <-  c("cbind(AD, BP) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      } else {
        null_model_formula <- c("cbind(AD, BP) ~ SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        alt_model_formula <-  c("cbind(AD, BP) ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      }
      # fitting the model
      position_error <- c()
      p_res_list <- list()
      for (p in seq_along(positions)){
          position <- positions[p]
          cat(paste0("\nOut of ",length(positions), " positions tested the following number has started: ",p,"\n"))
          pos_input_data <- 
            input_data %>%
              filter(Pos == position)
          # omit positions which throw errors
          skip_to_next <- FALSE
          # learning the reduced model
          tryCatch(
            invisible(
              capture.output(
                suppressWarnings(
                beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                                      ~ 1, 
                                      data = pos_input_data)
                          )
                        )
                      ) 
          , error = function(e) { skip_to_next <<- TRUE}
          )
          # learn the full model      
          tryCatch(
            invisible(
              capture.output(
                suppressWarnings(
              beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                                          ~ 1,
                                      data = pos_input_data)
                )
                        )
                      ) 
          , error = function(e) { skip_to_next <<- TRUE}
          )
          if(skip_to_next == FALSE){
          # log warnings
          convergence_warning_m0 <- FALSE
          convergence_warning_m1 <- FALSE
          singular_hessian_m0 <- FALSE
          singular_hessian_m1 <- FALSE
          if(beta_bin_m0@code != 0){ convergence_warning_m0 <- TRUE }
          if(beta_bin_m1@code != 0){ convergence_warning_m1 <- TRUE }
          if(beta_bin_m0@singular.hessian != 0){ singular_hessian_m0  <- TRUE }
          if(beta_bin_m1@singular.hessian != 0){ singular_hessian_m1  <- TRUE }
          na_deviance_m0 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE) 
          na_deviance_m1 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE)
          # extract summary from full model result
          # res <- summary(beta_bin_m1)@Coef[2, ]
          # 2023_06_06: there's some updates here/new stuff here calls of the result 
          # object and hence the way to access it changed. 
          # https://search.r-project.org/CRAN/refmans/aod/html/summary.glimML-class.html
          # although at first one might think it's an S3 vs S4 type object issue
          # in fact it is a class issue glimML 
          # --> hack which solved it: use show() instead of summary()  
          res <- show(beta_bin_m1)
          # assemble result data frame
          tmp_df <- data.frame(Pos = position, 
                              n_samples = length(unique(pos_input_data$SUBJID)),
                              estimate = res@Coef[2, "Estimate"], 
                              std_error = res@Coef[2, "Std. Error"], 
                              m0_p_value = res@Coef[1, "Pr(> |z|)"],
                              m1_p_value = res@Coef[2, "Pr(> |z|)"],
                              m1_z_value = res@Coef[2, "z value"],
                              intercept = res@Coef[ 1, "Estimate"],
                              intercept_p_value = res@Coef[ 1, "Pr(> |z|)"],
                              LR_m0 = beta_bin_m0@logL,
                              LR_m1 = beta_bin_m1@logL,
                              likelihood_ratio = 2*(beta_bin_m1@logL - beta_bin_m0@logL),
                              df_m0 =  beta_bin_m0@df.residual,
                              df_m1 =  beta_bin_m1@df.residual,
                              convergence_warning_m0 = convergence_warning_m0,
                              convergence_warning_m1 = convergence_warning_m1,
                              singular_hessian_m0 = singular_hessian_m0,
                              singular_hessian_m1 =  singular_hessian_m1,
                              na_deviance_m0 = na_deviance_m0,
                              na_deviance_m1 = na_deviance_m1,
                              iterations_m0 = beta_bin_m0@iterations,
                              iterations_m1 = beta_bin_m1@iterations,
                              stringsAsFactors = FALSE)
          # store result in list 
          p_res_list[[p]] <- tmp_df
        } else {
        position_error <- c(position_error,position)
      }
    }
    if(length(p_res_list) > 0){
      # collapse to data frame
      t_res_df <- bind_rows(p_res_list) 
      # record positions with errors
      t_res_df$pos_with_errors <- position_error
      # do LR test for getting analytical p-values
      t_res_df["LR_pval"] <- pchisq(t_res_df[ ,"likelihood_ratio"], df = 1, lower.tail = FALSE, log.p = FALSE) 
      # add tissue
      t_res_df$tissue <- tissue
      # retrieve runtime
      end_time <- Sys.time()
      t_res_df$run_time_seconds <- round(abs(difftime(start_time, end_time, units = "secs",)),3)
      print(paste0("Run time for tissue ",tissue," was ",unique(t_res_df$run_time_seconds)," seconds."))
      # 
      # store results
      t_res_list[[t]] <- t_res_df
    }
  }
}
BB_A_donor_age_df <- bind_rows(t_res_list)
# do multiple testing correction
BB_A_donor_age_df$BB_A_FDR <- p.adjust(BB_A_donor_age_df$LR_pval, method = 'fdr')
BB_A_donor_age_df$BB_A_bonferroni <- p.adjust(BB_A_donor_age_df$LR_pval, method = 'bonferroni')


#+ save results to disc --------------------------------------------------------
write.table(BB_A_donor_age_df,paste0(g_ops$results_dir,"2024_07_23_gt_BB_A_donor_age_all_tissues_common_het_pos.tsv"))

