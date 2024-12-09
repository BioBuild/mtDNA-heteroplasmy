# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library(data.table)
metadata_dir <- "/cis_eQTL/"
permutations_dir <- "/cis_eQTL/per_tissue/"
permutations_pooled_dir <- "/cis_eQTL/"

LM_analytical_pvals_df <- read.table(paste0(metadata_dir,"2024_07_23_gt_cis_eQTL_LM_A_all_tissues_poly_A_peer_factors_regressed_out.tsv"))

# redo calculation of bonferroni and FDR - this time study wide
LM_analytical_pvals_df <- 
    LM_analytical_pvals_df %>%
        mutate(analytical_FDR = p.adjust(LR_pval, method = "fdr"),
               analytical_bonferroni = p.adjust(LR_pval, method = "bonferroni"))
# poly A genes as defined in Rackham et al. 2022, Nat. Rev. Genet.
non_poly_A_protein_coding_gene <- c("MT-ND6")
poly_A_genes_IDs <- g_ops$mt_regions %>% 
    filter(transcript_type == "protein_coding",
           !gene == non_poly_A_protein_coding_gene) %>% 
    pull(ID)
# annotate if poly A gene
LM_analytical_pvals_df <- LM_analytical_pvals_df %>%
    mutate(poly_A = case_when(feature_id %in% poly_A_genes_IDs ~ TRUE, 
                              TRUE ~ FALSE))
# now do the subsetting
LM_A_sig_poly_A_only_df <- 
    LM_analytical_pvals_df %>% 
        filter(analytical_bonferroni <= 0.05, poly_A == TRUE) %>%
        mutate(keep = TRUE) %>%
        dplyr::select(tissue,Pos,feature_id,keep)

BB_aod_analytical_pvals_df <- readRDS(paste0(metadata_dir,"2024_07_23_gt_cis_eQTL_BB_A_all_tissues_lm_bonferroni_pos_peer_factors_regressed_out.rds"))
BB_aod_analytical_pvals_df <- 
    left_join(BB_aod_analytical_pvals_df,
              LM_A_sig_poly_A_only_df) %>%
    filter(keep == TRUE)

#+ loop through tissues and pool permutations ----------------------------------
# only consider tissues for which BB_A has been run
tissues <- unique(BB_aod_analytical_pvals_df$tissue)
filename_BB_permutations <- paste0(permutations_pooled_dir,"2024_07_23_gt_BB_permutations_cis_eQTL_all_tissues.rds")
filename_tissues_no_perms <- paste0(metadata_dir,"2024_07_23_gt_list_tissues_for_wich_no_BB_permutations.txt")
if(!file.exists(filename_BB_permutations)){
    permutations_list <- list()
    tissues_no_permutations <- c()
    for(t in seq_along(tissues)){
        # set tissue
        tissue <- tissues[t]
        cat(paste0("\ntissue #",t,"/",length(tissues),": ",tissue," ------------ \n"))
        # pool permutation files (name is a bit misleading; it's 50 files in fact)
        #+ save results to disc --------------------------------------------------------
        tissue_parallel_jobs <- list.files(path = paste0(permutations_dir,tissue,"/"), pattern = paste0("2024_06_13_gt_",tissue,"_cis_eQTL_peer_factors_regressed_out_bb_permutation_n_200_permutations_nr_parallel_*"))
        # check if they are all there 
        if(length(tissue_parallel_jobs) > 0){
        print(paste0(length(tissue_parallel_jobs), " permutation result files are present for tissue: ",tissue))
        # init list for looping
        tissue_parallel_jobs_list <- list()
        for(p in seq_along(tissue_parallel_jobs)){
            cat(paste0("\njobfile #",p,"/",length(tissue_parallel_jobs)))
            # now faster
            jobfile_df <- fread(paste0(g_ops$permutations_dir,tissue,"/",tissue_parallel_jobs[p]))
            succ_perms <- max(jobfile_df$permutation)
            jobfile_df$permutation <- jobfile_df$permutation + permutations_done_before 
            permutations_done_before <- sum(permutations_done_before,succ_perms)
            # store in list 
            tissue_parallel_jobs_list[[p]] <- jobfile_df
        }
        tissue_permutations_df <- bind_rows(tissue_parallel_jobs_list)
        tissue_permutations_df$tissue <- tissue
        permutations_list[[t]] <- tissue_permutations_df
        } else {
        tissues_no_permutations <- c(tissues_no_permutations,tissue)
        }
    }
    permutations_all_t_df <- bind_rows(permutations_list)
    # save to disc
    saveRDS(permutations_all_t_df,filename_BB_permutations)
    write.table(tissues_no_permutations,filename_tissues_no_perms)
} else {
    BB_aod_permutation_pvals_df <- readRDS(filename_BB_permutations)
}

## check number of permutations per tissue
BB_aod_permutation_pvals_df %>%
    distinct(tissue,permutation) %>%
    count(tissue) %>%
    arrange(desc(n)) %>% 
    as.data.frame()

## use 10000 permutations, if not enough run more
thr_n_perms <- 10000
BB_aod_permutation_pvals_df <- 
    BB_aod_permutation_pvals_df %>%
        filter(permutation <= thr_n_perms)

#+ BB model calculate empirical_p_values ---------------------------------------
# add BB_aod_analytical_pvals to permutation df for calculating empirical pvals
BB_aod_permutation_pvals_df <- 
    left_join(BB_aod_permutation_pvals_df,
              BB_aod_analytical_pvals_df %>% dplyr::select(tissue,Pos,feature_id, analytical_p_value = LR_pval),
              by = c("tissue","Pos","feature_id")) %>% 
    distinct()
# caclulate empirical p-vals
BB_aod_empirical_pvals_df <- 
   BB_aod_permutation_pvals_df %>%
        # get empirical_p_value
        group_by(tissue,Pos,feature_id,analytical_p_value,n_permutations) %>%
          summarise(n_permutation_p_vals_smaller = sum(LR_pval <= analytical_p_value)) %>%
        ungroup() %>%
          mutate(empirical_p_value = n_permutation_p_vals_smaller/n_permutations) %>%
        distinct()
# add them to analytical p-value data frame
BB_aod_empirical_pvals_df <-
    left_join(BB_aod_analytical_pvals_df,
              BB_aod_empirical_pvals_df %>% 
                dplyr::select(tissue,Pos,feature_id,n_permutations,n_permutation_p_vals_smaller,empirical_p_value),
              by = c("tissue","Pos","feature_id")) %>%
      distinct()
# multiple testing correction will be done downstream
 

#+ save result to disc ---------------------------------------------------------
write.table(BB_aod_empirical_pvals_df,paste0(metadata_dir,"2024_07_23_gt_BB_E_cis_eQTL_all_tissues_peer_factors_regressed_out.tsv"))

