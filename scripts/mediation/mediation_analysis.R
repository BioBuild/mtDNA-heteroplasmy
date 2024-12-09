# author:   simon.wengert@helmholtz-munich.de
source("~/scripts/utils/global_settings.R")
library("tidyverse")
overlap_df <- read_csv("2024_06_16_overlap_donor_age_cis_eQTL_sig_hits__peer_factors_regressed_out.csv")
genotype_path <- "/genotype_files/split_by_tissue/"

#+ explore mediation by donor age effect analysis ------------------------------
# let's do this for all the overlaps
load_these_df <- overlap_df %>% 
  distinct(tissue,Pos,feature_id) %>% 
  as.data.frame()
tissues <- unique(load_these_df$tissue)
tissue_res_list <- list()
for(t in seq_along(tissues)){
  # select tissue 
  tissue <- tissues[t]
  # get positions per tissue where there is an overlapping effect
  positions <- unique(load_these_df[which(load_these_df$tissue == tissue), "Pos"])
  # load genotypes
  genotypes_df <- read_tsv(paste0(genotype_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"))
  # init looping list
  feature_res_list <- list()
  for(f in seq_along(features)){
    # get feature
    feature <- features[f]
    # load gene_expression values 
    ge_file_df <- read_tsv(paste0("/Users/simon.wengert/data/mtDNA_variants/metadata/gene_expression/",tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt")) %>% 
        filter(feature_id == feature) %>%
      pivot_longer(-"feature_id",names_to = "SUBJID", values_to = "logTPM")
    # merge into one
    model_input_df <- left_join(genotypes_df,ge_file_df, by = "SUBJID") %>% drop_na()
    # init looping list
    pos_res_list <- list()
    for(p in seq_along(positions)){
      # subset data
      position <- positions[p]
      pos_input_df <- model_input_df[which(model_input_df$Pos == position), ]
      
      # let's standarise the variables of interest: 
      # heteroplasmy, donor age and gene expression
      pos_input_df$sum_heteroplasmic_level <- scale(pos_input_df$sum_heteroplasmic_level)
      pos_input_df$AGE <- scale(pos_input_df$AGE)
      pos_input_df$logTPM <- scale(pos_input_df$logTPM)
      
      # do partial correlation analysis
      #
      # residuals edge 1:
      ## get LMs
      lm_het_on_ge <- lm(sum_heteroplasmic_level ~ logTPM + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      lm_age_on_ge <- lm(AGE ~ logTPM + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      ## get residuals
      residuals_het_on_ge <- residuals(lm_het_on_ge)
      residuals_age_on_ge <- residuals(lm_age_on_ge)
      cor_test_i <- cor.test(residuals_het_on_ge,residuals_age_on_ge,method = "pearson")
      cor_test_i_coef <- cor_test_i$estimate
      cor_test_i_pval <- cor_test_i$p.value
      # keep a record of residuals for plotting
      edge_1_df <- tibble(tissue = tissue,
                          Pos = position,
                          feature_id = feature,
                          residuals_x = residuals_het_on_ge,
                          residuals_y = residuals_age_on_ge,
                          x_type = "residuals_het_on_ge",
                          y_type = "residuals_age_on_ge",
                          p_val = cor_test_i_pval,
                          coef = cor_test_i_coef,
                          correlation_type = "pearson",
                          edge_number = "1")      
      # make plot 
      if(F){
      plot(residuals_het_on_ge,residuals_age_on_ge, font.lab = 2, cex.lab = 1.3,
          main = paste0("Person correlation:\ncoef =  ",as.character(cor_test_i_coef),
                       "\npval = ",as.character(cor_test_i_pval)))
      }
      
      # residuals edge 2:
      ## get LMs
      lm_het_on_age <- lm(sum_heteroplasmic_level ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      lm_ge_on_age <- lm(logTPM ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      ## get residuals
      residuals_het_on_age <- residuals(lm_het_on_age)
      residuals_ge_on_age <- residuals(lm_ge_on_age)
      # get correlations
      cor_test_ii <- cor.test(residuals_het_on_age,residuals_ge_on_age,method = "pearson")
      cor_test_ii_coef <- cor_test_ii$estimate
      cor_test_ii_pval <- cor_test_ii$p.value
      # keep a record of residuals for plotting
      edge_2_df <- tibble(tissue = tissue,
                          Pos = position,
                          feature_id = feature,
                          residuals_x = residuals_het_on_age,
                          residuals_y = residuals_ge_on_age,
                          x_type = "residuals_het_on_age",
                          y_type = "residuals_ge_on_age",
                          p_val = cor_test_ii_pval,
                          coef =  cor_test_ii_coef,
                          correlation_type = "pearson",
                          edge_number = "2") 
      # make plot
      if(F){
      plot(residuals_het_on_age,residuals_ge_on_het, font.lab = 2, cex.lab = 1.3,
           main = paste0("Person correlation:\ncoef =  ",as.character(cor_test_ii_coef),
                         "\npval = ",as.character(cor_test_ii_pval)))
      }
      
      # residuals edge 3:
      ## get LMs
      lm_age_on_het <-  lm(AGE ~ sum_heteroplasmic_level + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      lm_ge_on_het <- lm(logTPM ~ sum_heteroplasmic_level + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = pos_input_df)
      ## get residuals
      residuals_age_on_het <- residuals(lm_age_on_het)
      residuals_ge_on_het <- residuals(lm_ge_on_het)
      # get correlations
      cor_test_iii <- cor.test(residuals_age_on_het,residuals_ge_on_het,method = "pearson")
      cor_test_iii_coef <- cor_test_iii$estimate
      cor_test_iii_pval <- cor_test_iii$p.value
      # keep a record of residuals for plotting
      edge_3_df <- tibble(tissue = tissue,
                          Pos = position,
                          feature_id = feature,
                          residuals_x = residuals_age_on_het,
                          residuals_y = residuals_ge_on_het,
                          x_type = "residuals_age_on_het",
                          y_type = "residuals_ge_on_het",
                          p_val = cor_test_iii_pval,
                          coef = cor_test_iii_coef,
                          correlation_type = "pearson",
                          edge_number = "3") 
      # make plot
      if(F){
      plot(residuals_age_on_het,residuals_ge_on_het,font.lab = 2, cex.lab = 1.3,
           main = paste0("Person correlation:\ncoef =  ",as.character(cor_test_iii_coef),
                         "\npval = ",as.character(cor_test_iii_pval)))
      }
      
      # build result data frame
      ## residuals and cortest results for further analysis and plotting
      pos_res_df <- bind_rows(edge_1_df,edge_2_df,edge_3_df)
      ## previous
      if(F){
      ## summary of correlation test results
      pos_res_df <- tibble(tissue = tissue,
                           Pos = position,
                           feature_id = feature,
                           cor_test_i_coef = cor_test_i_coef,
                           cor_test_i_pval = cor_test_i_pval,
                           cor_test_ii_coef = cor_test_ii_coef,
                           cor_test_ii_pval = cor_test_ii_pval,
                           cor_test_iii_coef = cor_test_iii_coef,
                           cor_test_iii_pval = cor_test_iii_pval)
     }
    
      pos_res_list[[p]] <- pos_res_df
    }
    feature_res_df <- bind_rows(pos_res_list)
    feature_res_list[[f]] <- feature_res_df
  }
  tissue_res_df <- bind_rows(feature_res_list)
  tissue_res_list[[t]] <- tissue_res_df
}
res_df <- bind_rows(tissue_res_list)
res_df <- left_join(res_df,load_these_df %>% mutate(keep =TRUE)) %>% filter(keep == TRUE)
res_summary_df <- res_df %>% distinct(tissue,Pos,feature_id,edge_number,x_type,y_type,correlation_type,p_val,coef)
# do bonferrroni for multiple testing correction
res_summary_df$bonferroni <- p.adjust(res_summary_df$p_val,method = "bonferroni")
res_summary_df$is_edge <- ifelse(res_summary_df$bonferroni <= 0.05, TRUE, FALSE)
# annotate full results data frame containing all residuals with bonferonni 
res_df <- 
  left_join(res_df,
            res_summary_df %>%
              dplyr::select(tissue,Pos,feature_id,edge_number,bonferroni,is_edge))


#+ save results ----------------------------------------------------------------
write_csv(res_df,"2024_06_16_donor_age_cis_eQTL_mediation_analyis_all_residuals_peer_factors_regressed_out.csv")
write_csv(res_summary_df,"2024_06_16_donor_age_cis_eQTL_mediation_analyis_summary_peer_factors_regressed_out.csv")

