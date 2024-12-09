# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")

LM_A_res_df <- read.table("/cis_eQTL/2024_07_23_gt_cis_eQTL_LM_A_all_tissues_poly_A_peer_factors_regressed_out.tsv")
BB_E_res_df <- read.table("/cis_eQTL/2024_07_23_gt_BB_E_cis_eQTL_all_tissues_peer_factors_regressed_out.tsv") %>% as_tibble()

#+ adding poly_A gene annotation -----------------------------------------------
# poly A genes as defined in Rackham et al. 2022, Nat. Rev. Genet.
non_poly_A_protein_coding_gene <- c("MT-ND6")
poly_A_genes_df <- g_ops$mt_regions %>% 
  filter(transcript_type == "protein_coding",
         !gene == non_poly_A_protein_coding_gene) %>% 
  distinct(feature_id = ID,gene_name = gene,poly_A = TRUE)

#+ perform STUDY WIDE multiple testing correction ------------------------------
# for linear model results
LM_A_res_df <- 
  LM_A_res_df %>%
       mutate(LM_A_FDR = p.adjust(LR_pval, method = 'fdr'),
              LM_A_bonferroni = p.adjust(LR_pval, method = 'bonferroni')) 
# for beta binomial regression results 
BB_E_res_df <- 
  BB_E_res_df %>%
        mutate(BB_A_FDR = p.adjust(LR_pval, method = 'fdr'),
               BB_A_bonferroni = p.adjust(LR_pval, method = 'bonferroni'),
               BB_E_FDR = p.adjust(empirical_p_value, method = 'fdr'),
               BB_E_bonferroni = p.adjust(empirical_p_value, method = 'bonferroni'))

#+ add genomic annotation ------------------------------------------------------
# subset to poly_A genes
LM_A_res_df <- 
  left_join(LM_A_res_df,
            poly_A_genes_df %>% dplyr::select(-poly_A),
            by = ("feature_id")) 


#+ join the above into one data frame for plotting -----------------------------
gt_cis_eqtl_res_df <- 
  left_join(LM_A_res_df %>%
              dplyr::select(tissue,Pos,feature_id,gene_name,poly_A,n_samples,
                            LM_estimate = estimate,LM_std_error = lm_std_error,
                            LM_A = LR_pval,LM_A_FDR,LM_A_bonferroni,),
            BB_E_res_df %>%
              dplyr::select(tissue,Pos,feature_id,n_permutations,
                            BB_estimate = estimate, BB_std_error = std_error,
                            BB_A = LR_pval,BB_A_FDR,BB_A_bonferroni,
                            BB_E_pval = empirical_p_value, BB_E_FDR, BB_E_bonferroni),
            by = c("tissue", "Pos", "feature_id")) %>%
  distinct()


#+ save to disc ----------------------------------------------------------------
write_tsv(gt_cis_eqtl_res_df, "/cis_eQTL/2024_07_23_gt_cis_eqtl_peer_factors_regressed_out_res_df.tsv")

