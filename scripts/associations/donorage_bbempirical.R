# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")

metadata_dir <- "/genotype_files/split_by_tissue/"
results_dir <- "/donor_age/"
permutations_dir <- "/donor_age/per_tissue/"

#+ load m1A_G_rna_meth annotation ----------------------------------------------
mt_tRNA_modifications_df <-
  read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv") %>%
  filter(rna_modification %in% c("m1A","m1G"))

#+ load results donor age analytical p-values ----------------------------------
BB_aod_analytical_pvals_df <- read.table(paste0(results_dir,"2024_07_23_gt_BB_A_donor_age_all_tissues_common_het_pos.tsv"))

#+ loop through tissues and pool permutations ----------------------------------
tissues <- unique(BB_aod_analytical_pvals_df$tissue)
filename_BB_permutations <- paste0(g_ops$results_dir,"2024_07_23_gt_BB_permutations_donor_age_all_tissues.rds")
if(!file.exists(filename_BB_permutations)){
    permutations_list <- list()
    for(t in seq_along(tissues)){
        # set tissue
        tissue <- tissues[t]
        cat(paste0("\ntissue #",t,"/",length(tissues),": ",tissue))
        # pool permutation files (name is a bit misleading; it's 50 files in fact)
        tissue_parallel_jobs <- list.files(path = paste0(g_ops$permutations_dir,tissue,"/"), pattern = paste0("2024_04_28_gt_",tissue,"_pheno_bb_permutation_n_150_permutations_job_nr*"))
        # not all of them will have worked out of experience
        print(paste0(length(tissue_parallel_jobs), " permutation result files are present for tissue: ",tissue))
        # now read in successful permutations  
        # init list for looping
        tissue_parallel_jobs_list <- list()
        permutations_done_before <- 0
        for(p in seq_along(tissue_parallel_jobs)){
            # read in batch file
            jobfile_df <- read_tsv(paste0(g_ops$permutations_dir,tissue,"/",tissue_parallel_jobs[p]), col_types = cols())
            # adjust n_permutations based on the batch number (because for each
            # patch ran in parallel another 50 permutations would have been run already)
            succ_perms <- max(jobfile_df$permutation)
            jobfile_df$permutation <- jobfile_df$permutation + permutations_done_before 
            # adjust permuataions count for next time
            permutations_done_before <- sum(permutations_done_before,succ_perms)
            # store in list 
            tissue_parallel_jobs_list[[p]] <- jobfile_df
        }
        tissue_permutations_df <- bind_rows(tissue_parallel_jobs_list)
        tissue_permutations_df$tissue <- tissue
        permutations_list[[t]] <- tissue_permutations_df
    }
    permutations_all_t_df <- bind_rows(permutations_list)
    # save to disc
    saveRDS(permutations_all_t_df,filename_BB_permutations)
} else {
    BB_aod_permutation_pvals_df <- readRDS(filename_BB_permutations)
}

thr_n_perms <- 10000

# 1) load linear model analytical runs 
LM_analytical_pvals_df <- read_tsv(paste0(results_dir,"2024_07_23_gt_LM_A_donor_age_all_tissues_common_het_pos.tsv"))
# 2) re-do multiple testing correction on analytical p-values: 
#    since this has been done per tisssue before - this time study-wide!!!
LM_analytical_pvals_df <- 
  LM_analytical_pvals_df %>%
    mutate(analytical_FDR = p.adjust(analytical_p_value, method = 'fdr'),
           analytical_bonferroni = p.adjust(analytical_p_value, method = 'bonferroni'))
# 3) filter 
BB_aod_permutation_pvals_df <- 
  left_join(BB_aod_permutation_pvals_df,
            LM_analytical_pvals_df %>% 
              filter(analytical_bonferroni <= 0.05) %>% 
              mutate(keep = TRUE) %>% 
              dplyr::select(Pos,tissue,keep),
            by = c("tissue","Pos")) %>%
    filter(keep == TRUE) %>%
    dplyr::select(-keep) %>%
    group_by(tissue) %>%
      filter(permutation <= thr_n_perms)
# 4) limit BB analytical pvalues data frame to the same set of positions
BB_aod_analytical_pvals_df <- 
  left_join(BB_aod_analytical_pvals_df,
            LM_analytical_pvals_df %>% 
              filter(analytical_bonferroni <= 0.05) %>% 
              mutate(keep = TRUE) %>% 
              dplyr::select(Pos,tissue,keep),
            by = c("tissue","Pos")) %>%
    filter(keep == TRUE) %>%
    dplyr::select(-keep) 


#+ BB model calculate empirical_p_values ---------------------------------------
BB_aod_permutation_pvals_df <- 
    left_join(BB_aod_permutation_pvals_df,
              BB_aod_analytical_pvals_df %>% dplyr::select(tissue,Pos, analytical_p_value = LR_pval),
              by = c("tissue","Pos")) %>% 
    distinct()
# choose max permutations run succ per tissue
BB_aod_permutation_pvals_df <- 
    BB_aod_permutation_pvals_df %>% 
        group_by(tissue) %>% 
            mutate(n_permutations = max(permutation)) %>% 
        ungroup()
# caclulate empirical p-vals
BB_aod_empirical_pvals_df <- 
   BB_aod_permutation_pvals_df %>%
        # get empirical_p_value
        group_by(tissue,Pos,analytical_p_value,n_permutations) %>%
          summarise(n_permutation_p_vals_smaller = sum(LR_pval <= analytical_p_value)) %>%
        ungroup() %>%
          mutate(empirical_p_value = n_permutation_p_vals_smaller/n_permutations) %>%
        distinct() %>%
        # multiple testing correction STUDY WIDE!
        mutate(empirical_FDR = p.adjust(empirical_p_value, method = 'fdr'),
               empirical_bonferroni = p.adjust(empirical_p_value, method = 'bonferroni')) 
# add them to analytical p-value data frame
BB_aod_empirical_pvals_df <-
    left_join(BB_aod_analytical_pvals_df,
              BB_aod_empirical_pvals_df %>% dplyr::select(tissue,Pos,n_permutations,n_permutation_p_vals_smaller,empirical_p_value,empirical_FDR,empirical_bonferroni),
              by = c("tissue","Pos")) %>%
      distinct()


#+ save result to disc ---------------------------------------------------------
write.table(BB_aod_empirical_pvals_df,paste0(results_dir,"2024_07_23_gt_BB_E_donor_age_all_tissues_common_pos.tsv"))
