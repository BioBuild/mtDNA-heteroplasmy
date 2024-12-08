
# author: simon.wengert@helmholtz-muenchen.de
#

#+ state command args from parent script ---------------------------------------
print(commandArgs(trailingOnly = TRUE))


#+ define beta bin function ----------------------------------------------------
beta_bin_calc_empirical_pval <- function(model_input_data = model_input_data,
                                         type_of_input_data = type_of_input_data,
                                         phenotype = phenotype,
                                         disease_status = disease_status,
                                         positions = positions,
                                         n_permutations = n_permutations,
                                         log_tmp_files_dir = log_tmp_files_dir,
                                         mc_cores = mc_cores,
                                         n_parallel = n_parallel){

#' fit beta-binomial regression model using empirical p-values
#'
#' Test for relationship of mtDNA heteroplasmic levels with donor phenotypes in
#' large cohort sequening data (i.e. GTEx and CommonMind consortia). Significance
#' is assessed by calculating emprical p-values via shuffeling donor phenotype
#' labels for n_permutations. Note: run time can be several days! Empirical p-values
#' are estimated by (i) fitting the model per mtDNA positon to the non-modified 
#' input data. This 'real' posion p-values, in our case from LR testing, are (ii) 
#' then compared with the p-values obtained by shuffling the donor phenotype 
#' lables for N times (n_permutations). By breaking any relationship between 
#' explanatory and response variable, this way the null distribution of p-values
#' in our data is estimated. (iii) Empirical p-values are then estimated as follows:
#' empirical_p_value = (n_permuted_pvalues < real_position_pvalue)/n_permutations
#'
#' @param model_input_data long-format ('tidy') data frame with the a minimum
#' of required columns: SUBJID, AGE, disease_status (0/1 binarised), SEX (0/1 
#' binarised), Pos (mtDNA_position), AD ("N sucesses": sum_heteroplasmic_level 
#' * sum_coverage (per mtDNA position)), BP ("N failures": sum_coverage - AD).
#' @param type_of_input_data declare if non-modified input data or simulated 
#' data, for assessing model calibration, was used as input for modelling.
#' @param phenotype phenotype relationship to test for (AGE, BP or SCZ).
#' @param disease_status declare if donors from healthy controls or disease group
#' were used as input (Control, BP or SCZ).
#' @param positions integer vector of mtDNA postions which should be tested.
#' @param n_permutation number of permutations used for estimating empirical
#' p-value distribution (usually 1000).
#'
#' @author Simon Wengert (simon.wengert@helmholtz-muenchen.de)
#'

  # dependencies
  require(aod)
  require(MASS)
  require(parallel)
  # set seed for label shuffling when calc empirical p-values
  # set.seed(12682) ---> not doing this any more (we don't want to be all 
  # the permutations to be the same!!!)

  # print data set to console or log file
  cat("\nSTART -------------------------------------------------------------\n") 
  print(paste0("Type of input data set is: ",type_of_input_data))

  # select formulas according to phenotype chosen
  if(disease_status %in% c("Control","BP","SCZ") && phenotype == "AGE"){

    # GTEx only donor_age testing which is why only here we need to check if
    # the tissue investigated is present in both sexes  
    if(length(unique(model_input_data$SEX)) > 1){   
      null_model_formula <- c("cbind(AD, BP) ~ SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(AD, BP) ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }else{
      null_model_formula <- c("cbind(AD, BP) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(AD, BP) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }
    cat(paste0("\nTesting for relationship with: ",phenotype," \nin disease_status: ",
                disease_status,"\nnull model used: ",null_model_formula,
                "\nalternative model used: ",alt_model_formula,"\n")) 
  
  }else if(disease_status %in% c("BP","SCZ") && phenotype %in% c("BP","SCZ")){
  
    null_model_formula <- c("cbind(AD, BP) ~ SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    alt_model_formula <-  c("cbind(AD, BP) ~ disease_status + SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    cat(paste0("\nTesting for relationship with: ",phenotype,"\nnull model used: ",
                null_model_formula,"\nalternative model used: ",alt_model_formula,"\n"))
  
  }else{
  
    cat(paste0("\nphenotype '",phenotype,"' or disease_status '",disease_status,
    "' or combination of which invalid.","\nselect valid phenotype ('AGE', 'BP' or 'SCZ') disease_status ('Control', 'BP' or 'SCZ') combination."))
    stopifnot(disease_status %in% c("Control","BP","SCZ") && phenotype %in% c("AGE", "BP", "SCZ"))
  
  }

  cat(paste0("\nFitting beta binomial regression for ",length(unique(model_input_data$Pos))," mtDNA positions.\n"))

  # define function fitting beta bin regression per position:
  fit_beta_bin_per_pos <- function(a){

    # filter for postion of interest
    pos_input_data <- 
        model_input_data %>%
          filter(Pos == positions[a])

    # fit_model
    beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                                ~ 1, 
                                data = pos_input_data)
    
    beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                                ~ 1,
                                data = pos_input_data)
    
    # log warnings
    convergence_warning_m0 <- FALSE
    convergence_warning_m1 <- FALSE
    singular_hessian_m0 <- FALSE
    singular_hessian_m1 <- FALSE
    # @code contains info about convergence from within betabin() function
    # see source code for details: https://github.com/cran/aod/blob/master/R/betabin.R 
    if(beta_bin_m0@code != 0){ convergence_warning_m0 <- TRUE }
    if(beta_bin_m1@code != 0){ convergence_warning_m1 <- TRUE }
    if(beta_bin_m0@singular.hessian != 0){ singular_hessian_m0  <- TRUE }
    if(beta_bin_m1@singular.hessian != 0){ singular_hessian_m1  <- TRUE }
    if(FALSE){
    if(beta_bin_m0@code != 0 | beta_bin_m1@code != 0){
      cat(paste0("\nFollowing warinings were produced at mtDNA position: ",positions[a],
                 "\nwarning null model:", 
                 "\nPossible convergence problem. Optimization process code: ", beta_bin_m0@code, " (see ?optim).",
                 "\nwarning alternative model:", 
                 "\nPossible convergence problem. Optimization process code: ", beta_bin_m1@code, " (see ?optim).\n"))
      pos_warning <- TRUE
    }
    }

    # return results
    res <- summary(beta_bin_m1)@Coef[2, ]
    return(data.frame(Pos = positions[a], 
                      n_samples = length(unique(pos_input_data$SUBJID)),
                      estimate = res[ , "Estimate"], 
                      std_error = res[ , "Std. Error"], 
                      m0_p_value = summary(beta_bin_m0)@Coef[ 1, "Pr(> |z|)"],
                      m1_p_value = res[ , "Pr(> |z|)"],
                      m1_z_value = res[ , "z value"],
                      intercept = summary(beta_bin_m1)@Coef[ 1, "Estimate"],
                      intercept_p_value = summary(beta_bin_m1)@Coef[ 1, "Pr(> |z|)"],
                      likelihood_ratio = beta_bin_m0@dev - beta_bin_m1@dev, 
                      convergence_warning_m0 = convergence_warning_m0,
                      convergence_warning_m1 = convergence_warning_m1,
                      singular_hessian_m0 = singular_hessian_m0,
                      singular_hessian_m1 =  singular_hessian_m1,
                      iterations_m0 = beta_bin_m0@iterations,
                      iterations_m1 = beta_bin_m1@iterations,
                      stringsAsFactors = FALSE))
  }
  
  # define function calculating empirical p-values beta bin regression:
  shuffle_donor_age_beta_bin_per_pos <- function(a){

    # filter for position of interest
    pos_input_data <- 
        model_input_data %>%
          filter(Pos == positions[a])
    # shuffle donor labels
    if(phenotype == "AGE"){
      pos_input_data$AGE <- sample(pos_input_data$AGE, replace = FALSE)  
    }else if(disease_status %in% c("BP","SCZ") && phenotype %in% c("BP","SCZ")){
      pos_input_data$disease_status <- sample(pos_input_data$disease_status, replace = FALSE)
    }

    # fit_model
    beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                                ~ 1, 
                                data = pos_input_data)
    
    beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                                ~ 1,
                                data = pos_input_data)
    # log warnings
    convergence_warning_m0 <- FALSE
    convergence_warning_m1 <- FALSE
    singular_hessian_m0 <- FALSE
    singular_hessian_m1 <- FALSE
    # @code contains info about convergence from within betabin() function
    # see source code for details: https://github.com/cran/aod/blob/master/R/betabin.R 
    if(beta_bin_m0@code != 0){ convergence_warning_m0 <- TRUE }
    if(beta_bin_m1@code != 0){ convergence_warning_m1 <- TRUE }
    if(beta_bin_m0@singular.hessian != 0){ singular_hessian_m0  <- TRUE }
    if(beta_bin_m1@singular.hessian != 0){ singular_hessian_m1  <- TRUE }
    # include more verbose output in case needed in the future 
    res <- summary(beta_bin_m1)@Coef[2, ]
    # return results
    return(data.frame(Pos = positions[a], 
                      n_samples = length(unique(pos_input_data$SUBJID)),
                      estimate = res[ , "Estimate"], 
                      std_error = res[ , "Std. Error"], 
                      m0_p_value = summary(beta_bin_m0)@Coef[ 1, "Pr(> |z|)"],
                      m1_p_value = res[ , "Pr(> |z|)"],
                      m1_z_value = res[ , "z value"],
                      intercept = summary(beta_bin_m1)@Coef[ 1, "Estimate"],
                      intercept_p_value = summary(beta_bin_m1)@Coef[ 1, "Pr(> |z|)"],
                      likelihood_ratio = beta_bin_m0@dev - beta_bin_m1@dev, 
                      convergence_warning_m0 = convergence_warning_m0,
                      convergence_warning_m1 = convergence_warning_m1,
                      singular_hessian_m0 = singular_hessian_m0,
                      singular_hessian_m1 =  singular_hessian_m1,
                      iterations_m0 = beta_bin_m0@iterations,
                      iterations_m1 = beta_bin_m1@iterations,
                      stringsAsFactors = FALSE))
  }


  #+ first bit: get "real"/non-permuted LR p-values per position -----------------
  cat("\nFIT BETA BIN REGRESSION - NO PERMUTATAIONS\n") 
  # apply beta bin regression per position on non permuted data
  pos_res_df <- do.call(rbind, lapply(seq_along(positions), FUN = fit_beta_bin_per_pos))


  #+ save analytical p-value as intermediate result ------------------------------
  # create log dir if file doesn't exist yet


  #+ second bit: calculate empirical p-values per pos by permuting donor ages ----
  cat("\nCALCULATE EMPIRICAL P-VALUES\n")
  permutations <- 1:n_permutations
  pos_res_age_shuffeld_n_perms_df <- do.call(rbind, mclapply(seq_along(permutations), FUN = function(n){
    # fit beta bin regression on shuffeld donor ages
    tmp_res <- do.call(rbind, lapply(seq_along(positions), FUN = shuffle_donor_age_beta_bin_per_pos))
    tmp_res$permutation <- n
    # append tmp result to a log file after each iteration for safety backup
    fout <- paste0(log_tmp_files_dir,"iterations_empirical_p_val_",
            type_of_input_data,"_",phenotype,"_",disease_status,"_",
            n_permutations,"_",n_parallel,".tsv")
    if(n == 1){
      fwrite(tmp_res, fout, append = F, sep="\t")
    } else {
      fwrite(tmp_res, fout, append = T, sep="\t")
    }
    # print progress to screen/log-file
    print(paste0("empirical p-value permutation completed is ", permutations[n], " out of ", length(permutations)))
    # return results
    return(data.frame(Pos = tmp_res[ , "Pos"], 
                      n_samples = tmp_res[ , "n_samples"],
                      estimate = tmp_res[ , "estimate"],
                      m0_p_value = tmp_res[ , "m0_p_value"],
                      m1_p_value = tmp_res[ , "m1_p_value"],
                      m1_z_value = tmp_res[ , "m1_z_value"],
                      intercept = tmp_res[, "intercept"],
                      intercept_p_value = tmp_res[ , "intercept_p_value"],
                      likelihood_ratio = tmp_res[ , "likelihood_ratio"],
                      convergence_warning_m0 = tmp_res[ , "convergence_warning_m0"],
                      convergence_warning_m1 = tmp_res[ , "convergence_warning_m1"],
                      singular_hessian_m0 = tmp_res[ , "singular_hessian_m0"],
                      singular_hessian_m1 = tmp_res[ , "singular_hessian_m1"],
                      iterations_m0 = tmp_res[ , "iterations_m0"],
                      iterations_m1 = tmp_res[ , "iterations_m1"],
                      permutation = n,  
                      stringsAsFactors = FALSE))
    }
   , mc.cores = mc_cores)
  )

  # stop track time 
  # end_time <- Sys.time()
  # print(paste0("Run time ", difftime(start_time, end_time, units = "hours")," hours."))

  # perform likelihood ratio test
  pos_res_age_shuffeld_n_perms_df["LR_pval"] <- pchisq(pos_res_age_shuffeld_n_perms_df[ ,"likelihood_ratio"], df = 1, lower.tail = FALSE, log.p = FALSE)
  
  # initialise NA column for empirical p-values 
  pos_res_df <-  
    pos_res_df %>%
      mutate(empirical_p_value = as.character(NA))
  
  # calculate empirical p-values per position
  for(a in seq_along(positions)){
    # 1. get the "real" p-value at this position from non permuted model
    pos_real_p_val <- pos_res_df %>% filter(Pos == positions[a]) %>% pull(LR_pval) 
    # 2. get number of permuted p-values which are smaller than the real p-value
    n_permuted_p_val_passing <- sum(pos_res_age_shuffeld_n_perms_df[which(pos_res_age_shuffeld_n_perms_df$Pos == positions[a]), "LR_pval"] <= pos_real_p_val)
    # 3. calculate empirical p-value for this position by: 
    pos_empirical_p_value <- ifelse(n_permuted_p_val_passing < 1, paste0("smaller_than_1_over_",n_permutations), as.character(n_permuted_p_val_passing/n_permutations))  
    # 4. calculate empirical p-value for this position by:
    pos_res_df <- 
      pos_res_df %>%
        mutate(empirical_p_value = case_when(Pos == positions[a] ~ as.character(pos_empirical_p_value),
                                             TRUE ~ empirical_p_value))
  }                                               

  # multiple testing correction
  # (i)  avoide deviding by character for empirical p values "too small to be seen"
  # (ii) avoid calculations for mtDNA pos for which no LRs could be computed
  #      due to weigths == 0 (see labjournal.md 2022_03_30 for details)
  # filter out both
  idx <- 
    pos_res_df %>%
      add_rownames() %>%
      filter(!empirical_p_value == paste0("smaller_than_1_over_",n_permutations),
             !is.na(LR_pval)) %>%
      pull(rowname) %>%
      as.numeric()
  pos_res_df[idx, "empirical_FDR"] <-
    p.adjust(unlist(pos_res_df[idx, "empirical_p_value"]), method = 'fdr')

  # "pooled" empirical p-values 
  # calculate empirical p-values across all positions assuming they come from 
  # the same distribution since same set of donors = same sample 
  idx_nans <- which(is.na(pos_res_df$LR_pval) == TRUE)
  nan_pos <- pos_res_df[idx_nans, ]$Pos
  idx_nans_permuted <- which(pos_res_age_shuffeld_n_perms_df$Pos %in% nan_pos)
  pooled_permuted_p_values <- pos_res_age_shuffeld_n_perms_df[-idx_nans_permuted , ]$LR_pval
  empirical_p_values_pooled <- 
    pos_res_df[-idx_nans,] %>%
      group_by(Pos) %>%
        mutate(empirical_p_value_pooled = 
                case_when(sum(pooled_permuted_p_values < LR_pval) == 0 ~ paste0("smaller_than_1_over_",length(pooled_permuted_p_values)),
                              TRUE ~ as.character(sum(pooled_permuted_p_values < LR_pval)/length(pooled_permuted_p_values)))) %>%
      ungroup() %>%
      pull(empirical_p_value_pooled)
  pos_res_df[-idx_nans, "empirical_p_value_pooled"] <- empirical_p_values_pooled
  # correct for multiple testing
  pos_res_df[-idx_nans, "empirical_FDR_pooled"] <-
    p.adjust(unlist(pos_res_df[-idx_nans, "empirical_p_value_pooled"]), method = 'fdr')
  
  # return result as list. the idea here is to save computationally expensive 
  # result of permutaitons as well. just in case.
  res_list <- list()
  res_list[[1]] <- pos_res_df
  res_list[[2]] <- pos_res_age_shuffeld_n_perms_df
  cat("\nEND ---------------------------------------------------------------\n") 
  return(res_list)                        
}













#+ below is still under development --------------------------------------------
# devnotes can be found in ~/git/projects/mtDNA_variants/scripts/2_phenotype_relationships/donor_age/gtex_v8/README.md
# main idea of the function below is to run permutation p-values first and only
# in order to quite pragmatically debug what happens to the permutations not run.



#+ define beta bin function ot run analytical p-values only --------------------
beta_bin_calc_analytical_pval <- function(model_input_data = model_input_data,
                                          type_of_input_data = type_of_input_data,
                                          phenotype = phenotype,
                                          disease_status = disease_status,
                                          positions = positions){               
  # dependencies
  require(aod)
  require(MASS)
  # set seed 
  set.seed(12682)

  # print data set to console or log file
  cat("\nSTART -------------------------------------------------------------\n") 
  print(paste0("Type of input data set is: ",type_of_input_data))

  # select formulas according to phenotype chosen
  if(disease_status %in% c("Control","BP","SCZ") && phenotype == "AGE"){

    # GTEx only donor_age testing which is why only here we need to check if
    # the tissue investigated is present in both sexes  
    if(length(unique(model_input_data$SEX)) > 1){   
      null_model_formula <- c("cbind(AD, BP) ~ SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(AD, BP) ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }else{
      null_model_formula <- c("cbind(AD, BP) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
      alt_model_formula <-  c("cbind(AD, BP) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }
    cat(paste0("\nTesting for relationship with: ",phenotype," \nin disease_status: ",
                disease_status,"\nnull model used: ",null_model_formula,
                "\nalternative model used: ",alt_model_formula,"\n")) 
  
  }else if(disease_status %in% c("BP","SCZ") && phenotype %in% c("BP","SCZ")){
  
    null_model_formula <- c("cbind(AD, BP) ~ SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    alt_model_formula <-  c("cbind(AD, BP) ~ disease_status + SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    cat(paste0("\nTesting for relationship with: ",phenotype,"\nnull model used: ",
                null_model_formula,"\nalternative model used: ",alt_model_formula,"\n"))
  
  }else{
  
    cat(paste0("\nphenotype '",phenotype,"' or disease_status '",disease_status,
    "' or combination of which invalid.","\nselect valid phenotype ('AGE', 'BP' or 'SCZ') disease_status ('Control', 'BP' or 'SCZ') combination."))
    stopifnot(disease_status %in% c("Control","BP","SCZ") && phenotype %in% c("AGE", "BP", "SCZ"))
  
  }

  cat(paste0("\nFitting beta binomial regression for ",length(unique(model_input_data$Pos))," mtDNA positions.\n"))

  # define function fitting beta bin regression per position:
  fit_beta_bin_per_pos <- function(a){

    # filter for postion of interest
    pos_input_data <- 
        model_input_data %>%
          filter(Pos == positions[a])

    # fit_model
    beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                                ~ 1, 
                                data = pos_input_data)
    
    beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                                ~ 1,
                                data = pos_input_data)
    
    # log warnings
    convergence_warning_m0 <- FALSE
    convergence_warning_m1 <- FALSE
    singular_hessian_m0 <- FALSE
    singular_hessian_m1 <- FALSE
    # @code contains info about convergence from within betabin() function
    # see source code for details: https://github.com/cran/aod/blob/master/R/betabin.R 
    if(beta_bin_m0@code != 0){ convergence_warning_m0 <- TRUE }
    if(beta_bin_m1@code != 0){ convergence_warning_m1 <- TRUE }
    if(beta_bin_m0@singular.hessian != 0){ singular_hessian_m0  <- TRUE }
    if(beta_bin_m1@singular.hessian != 0){ singular_hessian_m1  <- TRUE }
    # 2023_01_28: also log if deviance is NaN or not
    na_deviance_m0 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE) 
    na_deviance_m1 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE)
    if(FALSE){
    if(beta_bin_m0@code != 0 | beta_bin_m1@code != 0){
      cat(paste0("\nFollowing warinings were produced at mtDNA position: ",positions[a],
                 "\nwarning null model:", 
                 "\nPossible convergence problem. Optimization process code: ", beta_bin_m0@code, " (see ?optim).",
                 "\nwarning alternative model:", 
                 "\nPossible convergence problem. Optimization process code: ", beta_bin_m1@code, " (see ?optim).\n"))
      pos_warning <- TRUE
    }
    }

    # return results
    res <- summary(beta_bin_m1)@Coef[2, ]
    return(data.frame(Pos = positions[a], 
                      n_samples = length(unique(pos_input_data$SUBJID)),
                      estimate = res[ , "Estimate"], 
                      std_error = res[ , "Std. Error"], 
                      m0_p_value = summary(beta_bin_m0)@Coef[ 1, "Pr(> |z|)"],
                      m1_p_value = res[ , "Pr(> |z|)"],
                      m1_z_value = res[ , "z value"],
                      intercept = summary(beta_bin_m1)@Coef[ 1, "Estimate"],
                      intercept_p_value = summary(beta_bin_m1)@Coef[ 1, "Pr(> |z|)"],
                      LR_m0 = beta_bin_m0@logL,
                      LR_m1 = beta_bin_m1@logL, 
                      likelihood_ratio = beta_bin_m0@dev - beta_bin_m1@dev, 
                      convergence_warning_m0 = convergence_warning_m0,
                      convergence_warning_m1 = convergence_warning_m1,
                      singular_hessian_m0 = singular_hessian_m0,
                      singular_hessian_m1 =  singular_hessian_m1,
                      na_deviance_m0 = na_deviance_m0,
                      na_deviance_m1 = na_deviance_m1,
                      iterations_m0 = beta_bin_m0@iterations,
                      iterations_m1 = beta_bin_m1@iterations,
                      LR_pval_direct = anova(beta_bin_m0,beta_bin_m1,test = "LRT")@anova.table[ 2, 'P(> Chi2)'],
                      stringsAsFactors = FALSE))
  }
  
  #+ calculate analytical p-values per position --------------------------------
  cat("\nFIT BETA BIN REGRESSION - NO PERMUTATAIONS\n") 
  # apply beta bin regression per position on non permuted data
  pos_res_df <- do.call(rbind, lapply(seq_along(positions), FUN = fit_beta_bin_per_pos))
  # perform likelihood ratio test
  pos_res_df["LR_pval"] <- pchisq(pos_res_df[ ,"likelihood_ratio"], df = 1, lower.tail = FALSE, log.p = FALSE)
  # multiple testing correction (not done here any more -- see minimal postprocessing script)
  # pos_res_df["LR_FDR"] <- p.adjust(pos_res_df[ ,"LR_pval"], method = 'fdr')

  cat("\nEND ---------------------------------------------------------------\n") 
  return(pos_res_df)                        
}






#+ define beta bin function permutation p-values only --------------------------

beta_bin_permutation_pval <- function(model_input_data = model_input_data,
                                     type_of_input_data = type_of_input_data,
                                     phenotype = phenotype,
                                     disease_status = disease_status,
                                     positions = positions,
                                     n_permutations = n_permutations,
                                     mc_cores = mc_cores,
                                     n_parallel = n_parallel){


#+ dependencies ----------------------------------------------------------------
require(aod)
require(MASS)
require(parallel)


#+ print data set to console or log file ---------------------------------------
cat("\nSTART ---------------------------------------------------------------\n") 
print(paste0("Type of input data set is: ",type_of_input_data))

# select formulas according to phenotype chosen
if(disease_status %in% c("Control","BP","SCZ") && phenotype == "AGE"){

  # GTEx only donor_age testing which is why only here we need to check if
  # the tissue investigated is present in both sexes  
  if(length(unique(model_input_data$SEX)) > 1){   
    null_model_formula <- c("cbind(AD, BP) ~ SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    alt_model_formula <-  c("cbind(AD, BP) ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
  }else{
    null_model_formula <- c("cbind(AD, BP) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    alt_model_formula <-  c("cbind(AD, BP) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
  }
  cat(paste0("\nTesting for relationship with: ",phenotype," \nin disease_status: ",
              disease_status,"\nnull model used: ",null_model_formula,
              "\nalternative model used: ",alt_model_formula,"\n")) 

}else if(disease_status %in% c("BP","SCZ") && phenotype %in% c("BP","SCZ")){

  null_model_formula <- c("cbind(AD, BP) ~ SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
  alt_model_formula <-  c("cbind(AD, BP) ~ disease_status + SEX + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
  cat(paste0("\nTesting for relationship with: ",phenotype,"\nnull model used: ",
              null_model_formula,"\nalternative model used: ",alt_model_formula,"\n"))

}else{

  cat(paste0("\nphenotype '",phenotype,"' or disease_status '",disease_status,
  "' or combination of which invalid.","\nselect valid phenotype ('AGE', 'BP' or 'SCZ') disease_status ('Control', 'BP' or 'SCZ') combination."))
  stopifnot(disease_status %in% c("Control","BP","SCZ") && phenotype %in% c("AGE", "BP", "SCZ"))

}

#+ define funtion for shuffeling donor age labels and fitting the model --------
# define function calculating empirical p-values beta bin regression
# (will be used in function below):
shuffle_donor_age_beta_bin_per_pos <- function(a){

# filter for position of interest
pos_input_data <- 
    model_input_data %>%
      filter(Pos == positions[a])
# shuffle donor labels
if(phenotype == "AGE"){
  pos_input_data$AGE <- sample(pos_input_data$AGE, replace = FALSE)  
}else if(disease_status %in% c("BP","SCZ") && phenotype %in% c("BP","SCZ")){
  pos_input_data$disease_status <- sample(pos_input_data$disease_status, replace = FALSE)
}

# fit_model
beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                            ~ 1, 
                            data = pos_input_data)

beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                            ~ 1,
                            data = pos_input_data)
# log warnings
convergence_warning_m0 <- FALSE
convergence_warning_m1 <- FALSE
singular_hessian_m0 <- FALSE
singular_hessian_m1 <- FALSE

# @code contains info about convergence from within betabin() function
# see source code for details: https://github.com/cran/aod/blob/master/R/betabin.R 
if(beta_bin_m0@code != 0){ convergence_warning_m0 <- TRUE }
if(beta_bin_m1@code != 0){ convergence_warning_m1 <- TRUE }
if(beta_bin_m0@singular.hessian != 0){ singular_hessian_m0  <- TRUE }
if(beta_bin_m1@singular.hessian != 0){ singular_hessian_m1  <- TRUE }
# 2023_01_28: also log if deviance is NaN or not
na_deviance_m0 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE) 
na_deviance_m1 <- case_when(beta_bin_m0@dev == "NaN" ~ TRUE, TRUE ~ FALSE)
# include more verbose output in case needed in the future 
res <- summary(beta_bin_m1)@Coef[2, ]
# return results
return(data.frame(Pos = positions[a], 
                  n_samples = length(unique(pos_input_data$SUBJID)),
                  estimate = res[ , "Estimate"], 
                  std_error = res[ , "Std. Error"], 
                  m0_p_value = summary(beta_bin_m0)@Coef[ 1, "Pr(> |z|)"],
                  m1_p_value = res[ , "Pr(> |z|)"],
                  m1_z_value = res[ , "z value"],
                  intercept = summary(beta_bin_m1)@Coef[ 1, "Estimate"],
                  intercept_p_value = summary(beta_bin_m1)@Coef[ 1, "Pr(> |z|)"],
                  LR_m0 = beta_bin_m0@logL,
                  LR_m1 = beta_bin_m1@logL, 
                  likelihood_ratio = beta_bin_m0@dev - beta_bin_m1@dev, 
                  # 2023_01_28: after having a (super) large number of NaNs for 
                  # the deviances together with a warning mesage that at least
                  # one line had weights 0, I decided to directly perform the
                  # LRT in order to at least get a p-value out for further
                  # discussions - this is not ideal and has to be clarified still. 
                  LR_pval_direct = anova(beta_bin_m0,beta_bin_m1,test = "LRT")@anova.table[ 2, 'P(> Chi2)'],
                  convergence_warning_m0 = convergence_warning_m0,
                  convergence_warning_m1 = convergence_warning_m1,
                  singular_hessian_m0 = singular_hessian_m0,
                  singular_hessian_m1 =  singular_hessian_m1,
                  na_deviance_m0 = na_deviance_m0,
                  na_deviance_m1 = na_deviance_m1,
                  iterations_m0 = beta_bin_m0@iterations,
                  iterations_m1 = beta_bin_m1@iterations,
                  stringsAsFactors = FALSE))
}


#+ fit permutation p-value -----------------------------------------------------
cat(paste0("\nFitting beta binomial regression for ",length(unique(model_input_data$Pos))," mtDNA positions.\n"))
cat("\nCALCULATE PERMUTATION P P-VALUES VALUES ONLY\n")
permutations <- 1:n_permutations

permutations_df <- do.call(rbind, mclapply(seq_along(permutations), FUN = function(n){
#permutations_df <- do.call(rbind, lapply(seq_along(permutations), FUN = function(n){
  # fit beta bin regression on shuffeld donor ages
  tmp_res <- do.call(rbind, lapply(seq_along(positions), FUN = shuffle_donor_age_beta_bin_per_pos))
  tmp_res$permutation <- n
  # append tmp result to a log file after each iteration for safety backup
  fout <- paste0(log_tmp_files_dir,"iterations_permutation_p_val_",
          type_of_input_data,"_",phenotype,"_",disease_status,"_",
          n_permutations,"_",n_parallel,".tsv")
  if(n == 1){
    fwrite(tmp_res, fout, append = F, sep="\t")
  } else {
    fwrite(tmp_res, fout, append = T, sep="\t")
  }
  # print progress to screen/log-file
  print(paste0("permutation completed is ", permutations[n], " out of ", length(permutations)))
  # return results
  return(data.frame(Pos = tmp_res[ , "Pos"], 
                    n_samples = tmp_res[ , "n_samples"],
                    estimate = tmp_res[ , "estimate"],
                    m0_p_value = tmp_res[ , "m0_p_value"],
                    m1_p_value = tmp_res[ , "m1_p_value"],
                    m1_z_value = tmp_res[ , "m1_z_value"],
                    intercept = tmp_res[, "intercept"],
                    intercept_p_value = tmp_res[ , "intercept_p_value"],
                    LR_m0 = tmp_res[ , "LR_m0"],
                    LR_m1 = tmp_res[ , "LR_m1"], 
                    likelihood_ratio = tmp_res[ , "likelihood_ratio"],
                    convergence_warning_m0 = tmp_res[ , "convergence_warning_m0"],
                    convergence_warning_m1 = tmp_res[ , "convergence_warning_m1"],
                    singular_hessian_m0 = tmp_res[ , "singular_hessian_m0"],
                    singular_hessian_m1 = tmp_res[ , "singular_hessian_m1"],
                    na_deviance_m0 = tmp_res[ , "na_deviance_m0"],
                    na_deviance_m1 = tmp_res[ , "na_deviance_m1"],
                    iterations_m0 = tmp_res[ , "iterations_m0"],
                    iterations_m1 = tmp_res[ , "iterations_m1"],
                    permutation = n,  
                    LR_pval_direct = tmp_res[ , "LR_pval_direct"],
                    stringsAsFactors = FALSE))
  }
 , mc.cores = mc_cores)
)
cat("\nPERFORM LR TEST\n")
# perform likelihood ratio test to calculate permutaion pvalues 
permutations_df$LR_pval <- pchisq(permutations_df$likelihood_ratio, 
                                  df = 1, lower.tail = FALSE, log.p = FALSE)                                    
cat("\nCOMPLETE ------------------------------------------------------------\n") 
# return result
return(permutations_df)
}
