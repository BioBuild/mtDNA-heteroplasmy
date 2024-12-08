# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
metadata_dir <- "/donor_age/"

#+ load annotations ------------------------------------------------------------
# RNA modified sites
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))
# genomic annotation
chrM_anno <- read_tsv("~/metadata/annotations/gencode_v35_annotation_chrM.gtf") %>% as_granges()

#+ load donor age association results ------------------------------------------
LM_A_res_df <- read_tsv(paste0(metadata_dir,"2024_07_23_gt_LM_A_donor_age_all_tissues_common_het_pos.tsv"))
BB_E_res_df <- read.table(paste0(metadata_dir,"2024_07_23_gt_BB_E_donor_age_all_tissues_common_pos.tsv")) %>% as_tibble()

#+ make data frame for plottings -----------------------------------------------
plot_pheno_res_df <-
  left_join(LM_A_res_df %>%
              dplyr::select(tissue, Pos, LM_n_samples = n_donors, LM_estimate = estimate, 
                            LM_std_error = std_error, LM_A_pval = analytical_p_value,
                            LM_A_FDR = analytical_FDR, LM_A_bonferonni = analytical_bonferroni),
            BB_E_res_df %>% 
              dplyr::select(tissue, Pos, BB_n_samples = n_samples, n_permutations, 
                            BB_estimate = estimate, BB_std_error = std_error,
                            BB_A_pval = LR_pval, BB_A_FDR, 
                            BB_A_bonferroni, BB_E_pval = empirical_p_value,
                            BB_E_FDR = empirical_FDR, BB_E_bonferroni = empirical_bonferroni),
            by = c("tissue","Pos"))


#+ add molecular process and genomic annotations -------------------------------
plot_pheno_res_df <- 
  plot_pheno_res_df %>%
    mutate(molecular_process = case_when(Pos %in% m1A_G_methylations$genomic_pos ~ "rna_modification",
                                         TRUE ~ "dna_mutation")) 
plot_pheno_res_df <-
  join_overlap_left(plot_pheno_res_df %>% 
                       mutate(seqnames = 'chrM',
                              start = Pos,
                              end = Pos) %>%
                       as_granges(),
                     chrM_anno %>%
                       plyranges::select(gene_name)) %>%
                      as_tibble() %>%
                      distinct() %>%
  mutate(gene_name = case_when(is.na(gene_name) ~ "non_coding",
                               TRUE ~ gene_name))


#+ add additional codings ------------------------------------------------------
# colour for transcript type
plot_pheno_res_df <- 
  left_join(plot_pheno_res_df,
            g_ops$mt_regions %>%
              dplyr::select(gene_name = gene, transcript_type))
# sginificance labels for highlighting in plots
# plot_pheno_res_df$stars <- ifelse(plot_pheno_res_df$LM_A_bonferonni <= 0.05, "*", " ")
#plot_pheno_res_df$stars <- ifelse(plot_pheno_res_df$BB_A_bonferroni <= 0.05, "**", plot_pheno_res_df$stars) 
plot_pheno_res_df$stars <- ifelse(plot_pheno_res_df$BB_A_bonferroni <= 0.05, "#", " ") 


#+ save result -----------------------------------------------------------------
write_tsv(plot_pheno_res_df,paste0(metadata_dir,"2024_07_23_gt_donor_age_all_tissues_annotated.tsv"))

