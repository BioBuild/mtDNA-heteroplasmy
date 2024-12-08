# author:  simon.wengert@helmholtz-muenchen.de
# purpose: run beta binomial regression for donor age association testing.
#          Only run slow and intesive BB models on heteroplasmies significant 
#          in LM_analytical run. 


#+ specs and option parsing ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
machine <- args[1]
date <- args[2]
tissue <- args[3]
phenotype <- args[4]
disease_status <- args[5]
N_donors_min <- args[6]
heteroplasmy_threshold <- args[7]
n_permutations <- args[8]
mc_cores <- args[9]
cohort <- args[10]
n_parallel <- args[11]
# the latter is required in order to work with sourced scripts like the function
# below. See also here: https://stackoverflow.com/questions/14525580/how-to-pass
# -command-line-arguments-when-calling-source-on-an-r-file-within-ano
commandArgs <- args
# save dir per tissue
tissue_output <- paste0("/donor_age/per_tissue/",tissue,"/")
ifelse(!dir.exists(file.path(tissue_output)),dir.create(file.path(tissue_output),recursive = TRUE),FALSE)


#+ source dependencies ---------------------------------------------------------
system("conda activate gtex_mt_variants")
# source global options
source("~/scripts/utils/global_settings.R")
library("data.table")
# this custom function fits the beta-bin regression model
source('~/scripts/associations/functions_bb_empirical.R')


#+ R log files -----------------------------------------------------------------
log_dir <- "/donor_age/sink_logs/"
# for debugging
log_file <- paste0(log_dir,date,'_',tissue,'_gt_beta_bin_emp_pval_in_',disease_status,
                   '_donors_relationship_with_',phenotype,'_on_original_data_n_permutations_',
                   n_permutations,'_n_parallel_',n_parallel,'_DEBUG_file.txt')
# for partial permutation results in case slurm interruption
log_tmp_files_dir <- "/donor_age/perms_tmp/"
sink(log_file, append = TRUE)


#+ state parameters chosen -----------------------------------------------------
cat(paste0("\ncomputation performed on: ",machine,"\n"))
cat(paste0("\ncohort is: ",cohort,"\n"))
cat(paste0("\ndate is: ",date,"\n"))
cat(paste0("\ntissue is: ",tissue,"\n"))
cat(paste0("\nphenotype is: ",phenotype,"\n"))
cat(paste0("\ndisease_status is: ",disease_status,"\n"))
cat(paste0("\nN_donors_min is: ",f_ops$site_n_donors_min,"\n"))
cat(paste0("\nheteroplasmy_threshold is: ",heteroplasmy_threshold,"\n"))


#+ load data for respective tissue ---------------------------------------------
input_dir <- "/genotype_files/split_by_tissue/"
input_data <- read_tsv(paste0(input_dir,tissue,"/2023_07_10_gt_",tissue,"_heteroplasmy_genotypes_long_format.tsv"))
type_of_input_data <- c("original_data")
cat(paste0("\ntotal sample size is : ",length(unique(input_data$SUBJID)),"\n"))


#+ filter for significant heteroplasmy results from LM_analytical -------------- 
# we don't need to filter out pos with zero variance and with fewer than 60 
# donors since we already did that before when fitting the LM 
lm_res_all_tissues_df <- read_tsv("/donor_age/2023_07_10_gt_LM_A_donor_age_all_tissues.tsv")
pos_idx <- which(lm_res_all_tissues_df$tissue == tissue & lm_res_all_tissues_df$analytical_FDR <= 0.05)
positions <- as.data.frame(lm_res_all_tissues_df)[pos_idx, "Pos"]

################################################################################
###                     apply filters testing stage                          ###
################################################################################

#+ drop na's for now doing complete case analysis ------------------------------
input_data <- input_data %>% drop_na()


#+ filter out positions 0 variance ---------------------------------------------
input_data  <- 
  input_data %>%
    filter(var_per_pos > f_ops$site_variance_larger_than)


#+ filter for minimum number of donors -----------------------------------------
min_donors_pos <- 
  input_data  %>%
      group_by(Pos) %>%
      count(Pos) %>% 
        filter(n >= f_ops$site_n_donors_min) %>%
      ungroup() %>%  
      pull(Pos)
positions <- positions[which(positions %in% min_donors_pos)]


#+ prevent fitting model if no positions pass filter criteria ------------------
if(length(positions) == 0){
  stop(paste0("No mtDNA positions are passing selected filters for tissue: ", tissue))
} else {  


#+ prep model input ------------------------------------------------------------
relevant_columns <- c("SUBJID", "tissue", "Pos", "AD", "BP", "SEX", "AGE", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5")
model_input_data <- input_data[which(input_data$Pos %in% positions), relevant_columns]


#+ BB_permutation p_values -----------------------------------------------------
BB_P_start_time <- Sys.time()
# fitting the model
BB_P_res_df <- 
    beta_bin_permutation_pval(model_input_data = model_input_data,
                              type_of_input_data = type_of_input_data,
                              positions = as.integer(positions),
                              phenotype = phenotype,
                              n_permutations = n_permutations,
                              disease_status = disease_status,
                              mc_cores = mc_cores,
                              n_parallel = n_parallel,
                              log_tmp_files_dir = log_tmp_files_dir)
# calc runtime
BB_P_end_time <- Sys.time()
BB_P_res_df$run_time_hours <- abs(round(difftime(BB_P_start_time, BB_P_end_time, units = "hours"), digits = 5))
cat(paste0("\nRun time BB permutation p-values ", unique(BB_P_res_df$run_time_hours)," hours.\n"))
# save to disc
write_tsv(BB_P_res_df,paste0(tissue_output,"2023_07_11_gt_bb_permutation_n_",n_permutations,"_permutations_nr_parallel_",n_parallel,"_donor_age_",tissue,".tsv"))
}


#+ session_info ----------------------------------------------------------------
cat("\nSESSION INFO --------------------------------------------------------\n") 
sessionInfo()
sink()
