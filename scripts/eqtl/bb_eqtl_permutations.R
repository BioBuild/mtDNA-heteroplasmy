# author:  simon.wengert@helmholtz-munich.de

#+ option parsing --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
machine <- args[1]
cohort <- args[2]
date <- args[3]
tissue <- args[4]
n_permutations <- args[5]
mc_cores <- args[6]
n_parallel <- args[7]
signal_type <- args[8]
selection_mode <-args[9] 
commandArgs <- args

#+ set paths and hard coded params ---------------------------------------------
genotype_dir <- "/genotype_files/split_by_tissue/"
ge_dir <- "/gene_expression/"
tissue_output <- paste0("/cis_eQTL/per_tissue/",tissue,"/")
ifelse(!dir.exists(file.path(tissue_output)),dir.create(file.path(tissue_output),recursive = TRUE),FALSE)

#+ source dependencies ---------------------------------------------------------
source("~/scripts/utils/global_settings.R")
library("data.table")
source('~/scripts/associations/functions_bb_empirical.R')

#+ get filter options ----------------------------------------------------------
# doing this manually for now since sourcing global options gave conflicts on 
# clusters for some reason (2024-04-25).
all_filters_df <- read_csv("~/metadata/utils/all_filters.csv")
# make filter options list for direct access in pipeline
f_ops <- as.list(all_filters_df$value)
names(f_ops) <- all_filters_df$name
f_ops$sites_seq_artefacts <- as.integer(unlist(strsplit(f_ops$sites_seq_artefacts, " ")))
f_ops$site_signal_type_exploration <- unlist(strsplit(f_ops$site_signal_type_exploration, split = " "))
f_ops$variant_type <- as.integer(f_ops$variant_type)
f_ops$top_wmh <- as.double(f_ops$top_wmh)
f_ops$min_cov_fwd <- as.integer(f_ops$min_cov_fwd)
f_ops$min_cov_rev <- as.integer(f_ops$min_cov_rev)
f_ops$delta_frac <- as.double(f_ops$delta_frac)
f_ops$fraction_sum_of_levels_per_pos <- as.double(f_ops$fraction_sum_of_levels_per_pos)
f_ops$site_variance_larger_than <- as.integer(f_ops$site_variance_larger_than)
f_ops$het_n_donors_min <- as.integer(f_ops$het_n_donors_min)
f_ops$mad_scaled_median <- as.double(f_ops$mad_scaled_median)
f_ops$tissue_n_donors_min <- as.integer(f_ops$tissue_n_donors_min)
f_ops$donor_fraction_missing_sites <- as.double(f_ops$donor_fraction_missing_sites)
f_ops$sample_fraction_missing_sites <- as.double(f_ops$sample_fraction_missing_sites)
f_ops$site_n_donors_min <- as.integer(f_ops$site_n_donors_min)

#+ R console logging -----------------------------------------------------------
log_dir <- paste0("/cis_eQTL/sink_logs/",tissue,"/")
ifelse(!dir.exists(file.path(log_dir)),dir.create(file.path(log_dir),recursive = TRUE),FALSE)
log_file <- paste0(log_dir,date,'_',tissue,'_gt_BB_permutations_cis_eQTL_n_permutations_',
                   n_permutations,'_n_parallel_',n_parallel,'_DEBUG_file.txt')
# start logging
sink(log_file, append = TRUE)


#+ state parameters chosen -----------------------------------------------------
cat(paste0("\ncomputation performed on: ",machine,"\n"))
cat(paste0("\ncohort is: ",cohort,"\n"))
cat(paste0("\ndate is: ",date,"\n"))
cat(paste0("\ntissue is: ",tissue,"\n"))
cat(paste0("\n# donors min is: ",f_ops$site_n_donors_min,"\n"))
cat(paste0("\n# permutation is: ",n_permutations,"\n"))
cat(paste0("\nmc_cores is: ",mc_cores,"\n"))
cat(paste0("\nn_parallel is: ",n_parallel,"\n"))


#+ load input data per tissue --------------------------------------------------
# loading heteroplasmy genotypes
t_genotypes_df <- 
    read_tsv(paste0(genotype_dir,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv")) %>%
        dplyr::select(snp_id = feature_id,everything())
# loading gene expression files for poly_A genes
t_ge_file <- read_tsv(paste0(ge_dir,tissue,"_mtgenes_poly_A.txt"))
# make long to wide
t_ge_file <- melt(setDT(t_ge_file), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")

################################################################################
###   load LM_A QTL results and subset to study-wide LM_A bonferroni hits    ###
################################################################################
if(signal_type == "exploration"){
  # load linear model result
  #LM_A_all_tissues_df <- read.table("/cis_eQTL/2024_04_28_gt_cis_eQTL_LM_A_all_tissues_poly_A.tsv")
  LM_A_all_tissues_df <- read.table("/cis_eQTL/2024_06_13_gt_cis_eQTL_LM_A_all_tissues_poly_A_peer_factors_regressed_out.tsv")
  # recalculate bonferroni: tissue-wide before now we want to go with study wide
  LM_A_all_tissues_df <- LM_A_all_tissues_df %>% mutate(analytical_bonferroni = p.adjust(LR_pval, method = 'bonferroni'))
  # decide which set of position-feature pairs to calculate permutations for
  if(selection_mode == "lm_bonferroni_passing"){
    # pull out indices of significant positions - feature pairs
    test_pairs <- LM_A_all_tissues_df[which(LM_A_all_tissues_df$tissue == tissue & LM_A_all_tissues_df$analytical_bonferroni <= 0.05), c("Pos","feature_id")]
  } else if (selection_mode == "common_pos"){
    # keep all the positions that passed LM testing filters
    test_pairs <- LM_A_all_tissues_df[which(LM_A_all_tissues_df$tissue == tissue), c("Pos","feature_id")]
  }
} else if (signal_type == "confirmation"){
  BB_E_all_tissues_df <- read.table("/cis_eQTL/2023_09_12_gt_cis_eQTL_BB_A_all_tissues_confirmation.tsv")
  test_pairs <- BB_E_all_tissues_df[which(BB_E_all_tissues_df$tissue == tissue & BB_E_all_tissues_df$BB_A_bonferroni <= 0.05), ]
}

################################################################################
###          apply testing stage filters to heteroplasmy genotypes           ###
################################################################################
#+ filter for heteroplasmy signal_type you want to allow for -------------------
if(signal_type == "exploration"){
  t_genotypes_df <-
      t_genotypes_df %>% 
          filter(signal_type %in% f_ops$site_signal_type_exploration)
} else if (signal_type == "confirmation"){
  t_genotypes_df <-
      t_genotypes_df %>% 
          filter(signal_type %in% f_ops$site_signal_type_stringent)
} else {
  cat(paste0("\n singal type of heteroplasmies not specified. Please either set to 'exploration' or 'confirmation'\n"))
}


#+ filter out positions 0 variance ---------------------------------------------
t_genotypes_df <- 
  t_genotypes_df %>%
    filter(var_per_pos > f_ops$site_variance_larger_than)


#+ filter for minimum number of donors per site before association testing -----
min_donors_pos <- 
  t_genotypes_df  %>%
      group_by(Pos) %>%
      count(Pos) %>% 
        filter(n >= f_ops$site_n_donors_min) %>%
      ungroup() %>%  
      pull(Pos)
test_pairs <- test_pairs[which(test_pairs$Pos %in% min_donors_pos), ]


#+ prevent fitting model if no positions pass filter criteria ------------------
if(length(test_pairs$Pos) == 0){stop(paste0("No mtDNA positions are passing selected filters for tissue: ", tissue))}


################################################################################
###          filter molecular phenotype file for features of interest        ###
################################################################################
#+ build QTL input data and filter for use_positions ---------------------------
input_data <- left_join(t_genotypes_df,t_ge_file, by = "SUBJID")


#+ select relevant variables and covariates  -----------------------------------
relevant_columns <- c("SUBJID", "tissue", "snp_id","feature_id","logTPM","Pos", 
                      "AD", "BP", "SEX", "AGE", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5")


#+ apply joint filter for position feature pairs USING JOINS -------------------
qtl_input_data <- inner_join(test_pairs,input_data, by = c("Pos", "feature_id")) 
qtl_input_data <- qtl_input_data[ ,relevant_columns]


#+ state tissue and # samples used for qtl testing -----------------------------
cat(paste0("\n tissue is: ",tissue," # samples is : ",length(unique(t_genotypes_df$SUBJID)),"\n"))


################################################################################
###                   test SNP-feature pairs sequentially                    ###
################################################################################
# tracking run time  
BB_P_start_time <- Sys.time()
# let's keep it simple for now and do sequential testing per feauture.
# if that takes too much time, we can a think about parallelizing this as well.
features <- unique(qtl_input_data$feature_id) 
f_list <- list()
# do sequential testing
for(f in seq_along(features)){
    feature <- features[f]
    snps <- unique(qtl_input_data[qtl_input_data$feature_id == feature, "Pos"])
    # state feature to be tested
    cat(paste0("\n feature is: ",feature," (# ",f,"/",length(features),") ----------------------------------------- \n"))
    # subset data to feature oi
    feature_input_data <- qtl_input_data[which(qtl_input_data$feature_id == feature), ]
    # now let's fit BB_permutations for all the features per position
    int_df <- BB_permutations_cis_eQTL(model_input_data = feature_input_data,
                                       feature_id = feature,
                                       positions = snps,
                                       n_permutations = n_permutations,
                                       mc_cores = mc_cores,
                                       n_parallel = n_parallel)
    f_list[[f]] <- int_df
}
cat("\n BB RUNS COMPLETE ---------------------------------------------------\n") 


#+ collapse list to data frame and add run_time for permutations ---------------
BB_P_res_df <- do.call(rbind, f_list)
# add signal type to permutation runs
BB_P_res_df$signal_type <- signal_type
# track time
BB_P_end_time <- Sys.time()
BB_P_res_df$run_time_hours <- abs(round(difftime(BB_P_start_time, BB_P_end_time, units = "hours"), digits = 5))
cat(paste0("\nRun time BB_permutation p-values ", unique(BB_P_res_df$run_time_hours)," hours.\n"))


#+ adding LR test pvalues for permutations -------------------------------------
BB_P_res_df$LR_pval <- pchisq(BB_P_res_df$likelihood_ratio, df = 1, lower.tail = FALSE, log.p = FALSE)  
cat("\nCALCULATING P-VALUES FROM LL RATIOS COMPLETE ------------------------\n") 


#+ save results to disc --------------------------------------------------------
write_tsv(BB_P_res_df,paste0(tissue_output,date,"_gt_",tissue,"_cis_eQTL_peer_factors_regressed_out_bb_permutation_n_",n_permutations,"_permutations_nr_parallel_",n_parallel,"_",signal_type,"_",selection_mode,".tsv"))


#+ session_info ----------------------------------------------------------------
cat("\nSESSION INFO --------------------------------------------------------\n") 
sessionInfo()
sink()


