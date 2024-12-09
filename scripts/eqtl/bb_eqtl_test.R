# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")

metadata_dir <- "/genotype_files/split_by_tissue/"
gene_expression_dir <- "/gene_expression/"
results_dir <- "/cis_eQTL/"
single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
tissues <- list.dirs(metadata_dir, full.names = FALSE)[-1]
tissues <- tissues[-which(tissues %in% c("kidney_cortex"))]

#+ set signal_type to exploration or a confirmation run ------------------------
signal_type <- "exploration"

#+ load previous results for restricting model fit to significant hits ---------
if(signal_type == "exploration"){
    LM_A_all_tissues_df <- read.table("/cis_eQTL/2024_07_23_gt_cis_eQTL_LM_A_all_tissues_poly_A_peer_factors_regressed_out.tsv")
    res_previous_df <- LM_A_all_tissues_df[which(LM_A_all_tissues_df$analytical_bonferroni <= 0.05), ]
} else if (signal_type == "confirmation"){
    BB_E_all_tissues_df <- read_tsv("/cis_eQTL/2023_09_01_gt_cis_eqtl_res_df.tsv")
    res_previous_df <- BB_E_all_tissues_df[which(BB_E_all_tissues_df$BB_E_bonferroni <= 0.05), ]
}

#+ only use tissues which are present in previous result df --------------------
tissues <- tissues[which(tissues %in% res_previous_df$tissue)]

#+ loop through tissues and fit lm model per tissue and to previous hits -------
BB_A_start_time <- Sys.time()
t_res_list <- list()
for(t in seq_along(tissues)){
    tissue <- tissues[t]
    cat(paste0("\n# tissues is ",t,"/",length(tissues)," ----------------- \n"))
    # load genotypes
    t_genotypes_df <- 
        read_tsv(paste0(metadata_dir,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv")) %>% 
        select(snp_id = feature_id,everything())
    positions <- unique(t_genotypes_df$Pos)  
     # select for LM_A FDR hits only
    pos_idx <- which(res_previous_df$tissue == tissue)
    sig_positions <- unique(as.data.frame(res_previous_df)[pos_idx, "Pos"])
    # apply **testing stage filters** to genotypes
    # 1) filter for heteroplasmy signal_type you want to allow for 
    if(signal_type == "exploration"){
       t_genotypes_df <- t_genotypes_df %>% filter(signal_type %in% f_ops$site_signal_type_exploration)
    } else if (signal_type == "confirmation"){
       t_genotypes_df <- t_genotypes_df %>% filter(signal_type %in% f_ops$site_signal_type_stringent)
    }
    # 2) filter out positions 0 variance
    t_genotypes_df <- t_genotypes_df %>% filter(var_per_pos > f_ops$site_variance_larger_than)
    # 3) filter for minimum number of donors per site before association testing 
    min_donors_pos <- 
        t_genotypes_df  %>%
            group_by(Pos) %>%
                count(Pos) %>% 
                filter(n >= f_ops$site_n_donors_min) %>%
            ungroup() %>%  
        pull(Pos)
    use_positions <- sig_positions[which(sig_positions %in% min_donors_pos)]
    t_genotypes_df <- t_genotypes_df[which(t_genotypes_df$Pos %in% use_positions), ]
    # load gene expression for poly_A genes - regressed out peer factors
    t_ge_file <- read_tsv(paste0(gene_expression_dir,tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"))
    # make long to wide
    t_ge_file <- melt(setDT(t_ge_file), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
    # merge genotypes and gene expression for getting model input data
    input_data <- left_join(t_genotypes_df,t_ge_file, by = "SUBJID")
    # loop through each mtDNA position and test for
    # pairwise association with genes
    positions <- unique(input_data$Pos)
    features <- unique(input_data$feature_id)
    # set model formula + omit SEX as covariate in single sex tissues
    if(tissue %in% single_sex_tissues){   
        null_model_formula <- c("cbind(AD, BP) ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        alt_model_formula <-  c("cbind(AD, BP) ~ logTPM + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    } else {
        null_model_formula <- c("cbind(AD, BP) ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
        alt_model_formula <-  c("cbind(AD, BP) ~ logTPM + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5")
    }
    # fitting the model
    pos_res_list <- list()
    for(p in seq_along(positions)){
        position <- positions[p]
        cat(paste0("\n# positions is ",p,"/",length(positions),"\n"))
        pos_input_data <- input_data[which(input_data$Pos == position), ]
        feature_res_list <- list()
        for(f in seq_along(features)){
            feature <- features[f]
            feature_input_data <- pos_input_data[which(pos_input_data$feature_id == feature), ] 
            # fit BB_analytical for cis-eQTL
            beta_bin_m0 <- aod::betabin(as.formula(null_model_formula),      
                                        ~ 1, 
                                        data = feature_input_data)

            beta_bin_m1 <- aod::betabin(as.formula(alt_model_formula), 
                                        ~ 1,
                                        data = feature_input_data)
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
            # build result data frame
            res <- show(beta_bin_m1)
            # assemble result data frame
            tmp_df <- data.frame(tissue = unique(pos_input_data$tissue),
                                 Pos = position, 
                                 feature_id = feature,
                                 n_samples = length(unique(feature_input_data$SUBJID)),
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
            # store result
            feature_res_list[[f]] <- tmp_df
        }
        pos_res_list[[p]] <- bind_rows(feature_res_list)
    }
    t_res_df <- bind_rows(pos_res_list)
    # calculate analytical pvals
    t_res_df$LR_pval <- pchisq(t_res_df$likelihood_ratio, df = 1, lower.tail = FALSE, log.p = FALSE)  
    # return result
    t_res_list[[t]] <- t_res_df
}
# collapse list
BB_eqtl_analytical_df <- bind_rows(t_res_list)
# add signal type column such that we can identify what has been done later on
BB_eqtl_analytical_df$signal_type <- signal_type
# record total run time
BB_A_end_time <- Sys.time()
BB_eqtl_analytical_df$run_time_hours <- abs(round(difftime(BB_A_start_time,BB_A_end_time, units = "hours"), digits = 5))


#+ multiple testing correction -------------------------------------------------
BB_eqtl_analytical_df$BB_A_bonferroni <- p.adjust(BB_eqtl_analytical_df$LR_pval, method = "bonferroni")


#+ save to disc ----------------------------------------------------------------
saveRDS(BB_eqtl_analytical_df,paste0(g_ops$results_dir,"2024_07_23_gt_cis_eQTL_BB_A_all_tissues_lm_bonferroni_pos_peer_factors_regressed_out.rds"))

