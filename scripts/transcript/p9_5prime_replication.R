# author:   simon.wengert@helmholtz-munich.de
# pupose:   Ali et al 2020. tested correlations between inferred m1A/G methylation
#           levels and mRNA gene expression of the immediate 5' gene in cases
#           where there was an immediate 5' gene. This was their way of assessing
#           the (molecular) consequences of the varition in m1A/G methylation.
#           doi: https://doi.org/10.1038/s42003-020-0879-3

source("~/scripts/utils/global_settings.R")
library(lmtest)
library("kableExtra")

#+ define replication set -----------------------------------------------------
# position gene pairs tested by Ali et al. 2020. This is: p9 RNA methylation 
# sites (Pos) and the immediate 5'  upstream gene (gene_name) (see table 2)
pos_gene_df <-
  data.frame(
    Pos = c(1610,3238,4271,5520,8303,9999,10413,12146,15896),
    gene_name = c("MT-RNR1","MT-RNR2","MT-ND1","MT-ND2","MT-CO2","MT-CO3", 
                  "MT-ND3","MT-ND4","MT-CYB"),
    strand_name = c("H-Strand","H-Strand","H-Strand","H-Strand","H-Strand",
                    "H-Strand","H-Strand","H-Strand","H-Strand"),
    tissue = "whole_blood",
    beta_gtex_ali_et_al = c(0.09739377,2.511323,0.2212783,0.03422895,0.005807783,-0.1192088,0.05458773,0.1370204,0.2027844),
    p_value_gtex_ali_et_al = c(0.4027136,0.0000216,0.1716995,0.6816734,0.9465179,0.04399886,0.31342,0.05093543,0.2833017),
    sig_in_gtex_ali_et_al = c(FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)) %>%
  left_join(., g_ops$mt_regions %>% 
              distinct(feature_id = ID, gene_name = gene,poly_A)) %>%
  select(Pos,gene_name,feature_id,poly_A,everything())
# we only consider poly_A genes which means that we do not consider MT-RNR1 or MT-RNR2
pos_gene_df <- pos_gene_df %>% filter(poly_A == TRUE)

# show genes and positions we use in testing 
pos_gene_df %>%
  kbl() %>%
  kable_paper(bootstrap_options = "striped", full_width = F)

#+ save pos_gene_df for further analysis ---------------------------------------
write_tsv(pos_gene_df,"2024_11_11_positions_tested_in_Ali_et_al_2020.tsv")

#+ source aux files ------------------------------------------------------------
# gtex sample lookup file
gt_v8_lookup_df <- read_tsv("gt_v8_lookup_df.tsv",show_col_types = FALSE)
# donor gPCs
gt_gPCs <- read_tsv("~/git/mtDNA_variants/metadata/annotations/gt_v8_gPCs_europ.tsv")


#+ load read_counts heteroplasmy and gene expression data ----------------------
read_level_summaries_df <-
  read_tsv("~/data/mtDNA_variants/results/7_transcript_processing/tissue_summaries/2024-10-08_all_tissues_sample_summary_from_read_level_analysis_filtered.tsv",
           show_col_types = FALSE) %>%
  filter(Pos %in% pos_gene_df$Pos) 


#+ build model input data ------------------------------------------------------
model_input_df <- 
  read_level_summaries_df %>% 
  distinct(biospecimen_repository_sample_id,tissue,Pos,
           heteroplasmic_level_mtDNA_server,read_counts_heteroplasmic_level) %>%
  left_join(.,
            gt_v8_lookup_df %>% 
              distinct(biospecimen_repository_sample_id,tissue,SUBJID,AGE,SEX),
            by = c("biospecimen_repository_sample_id", "tissue")) %>%
  left_join(.,
            gt_gPCs,
            by = "SUBJID")


#+ run association test for all tissues ----------------------------------------
# init result object
res_list <- list()
# loop through tissues
tissues <- unique(model_input_df$tissue)
tissues <- tissues[-which(tissues == "kidney_cortex")]
for(t in seq_along(tissues)){
  # get tissue t
  tissue <- tissues[t]
  # tell us where we are
  cat(paste0("\n tissue is #",t,"/",length(tissues)))
  # subset to tissue
  t_df <- model_input_df[which(model_input_df$tissue == tissue), ]
  # add in tissue gene expression file
  t_gene_expression <- read_tsv(paste0("/gene_expression/",tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"), show_col_types = FALSE)
  t_gene_expression <- melt(setDT(t_gene_expression), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  # join tissue gene expression file to model input data
  t_df <-left_join(t_df,t_gene_expression, by = "SUBJID") 
  # only consider the 7 position feature pairs where poly_A mRNA is immediately
  # upstream of p9 methylation site
  t_df <- t_df %>% filter(paste(Pos, feature_id) %in% paste(pos_gene_df$Pos, pos_gene_df$feature_id))
  # loop through positions
  positions <- unique(t_df$Pos)
  p_res_list <- list()
  for(p in seq_along(positions)){
    # get position p
    position <- positions[p]
    # subset to position
    p_df <- t_df[which(t_df$Pos == position), ]
    # loop through features
    features <- unique(p_df$feature_id)
    f_res_list <- list()
    for(f in seq_along(features)){
      # get feature f
      feature <- features[f]
      # subset to feature
      f_df <- p_df[which(p_df$feature_id == feature), ]
      # run association test
      single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
      if(unique(p_df$tissue) %in% single_sex_tissues){   
        lm_m0 <- lm(read_counts_heteroplasmic_level ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
        lm_m1 <- lm(read_counts_heteroplasmic_level ~ logTPM + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
      } else {
        lm_m0 <- lm(read_counts_heteroplasmic_level ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
        lm_m1 <- lm(read_counts_heteroplasmic_level ~ logTPM + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
      }
      # get results and do statistical test
      res_m1 <- summary(lm_m1)
      LR_test <- lrtest(lm_m1,lm_m0) 
      # build res table
      tmp_df <- data.frame(tissue = tissue,
                           Pos = position,
                           feature_id = feature,
                           n_samples_read_counts = length(unique(p_df$SUBJID)),
                           property_tested = "read_level_counts",
                           estimate_read_counts =  lm_m1$coefficients["logTPM"],
                           lm_std_error_read_counts = res_m1$coefficients["logTPM","Std. Error"],
                           LR_pval_read_counts =  LR_test[2,  'Pr(>Chisq)'])
      
      # return result  
      f_res_list[[p]] <- tmp_df
    }
    # return result 
    f_res_df <- bind_rows(f_res_list)
    p_res_list[[p]] <- f_res_df
  }
  # return result
  p_res_df <- bind_rows(p_res_list)
  res_list[[t]] <- p_res_df
}
# collapse list to data frame
res_df <- bind_rows(res_list)
# multiple testing correction: 
## study wide
res_df$bonferroni_read_counts <- p.adjust(res_df$LR_pval_read_counts,method = "bonferroni") 
res_df$BH_read_counts <- p.adjust(res_df$LR_pval_read_counts,method = "BH") 
## within tissue (as this seems closest to what Ali et al. had done.)
res_df <- res_df %>%
  group_by(tissue) %>%
  mutate(t_bonferroni_read_counts = p.adjust(LR_pval_read_counts, method = "bonferroni"),
         t_BH_read_counts = p.adjust(LR_pval_read_counts, method = "BH")) %>%
  ungroup()
## lets assign significance based on tissue wide bonferroni
res_df <- res_df %>% 
  mutate(sig_read_counts = case_when(t_bonferroni_read_counts <= 0.05 ~ TRUE,TRUE ~ FALSE))

#+ compare our results with ali et al ------------------------------------------
compare_df <- left_join(pos_gene_df,res_df) %>% mutate(Pos = as.factor(Pos)) %>%
  select(tissue,Pos,n_samples_read_counts,gene_name,feature_id,poly_A,beta_gtex_ali_et_al,
         estimate_read_counts,lm_std_error_read_counts,p_value_gtex_ali_et_al,
         LR_pval_read_counts,t_bonferroni_read_counts,sig_read_counts,
         sig_in_gtex_ali_et_al)

# table
compare_df %>%
  kbl() %>%
  kable_paper(bootstrap_options = "striped", full_width = F)

# save whole blood replication overview
write_tsv(compare_df,"~/data/mtDNA_variants/results/7_transcript_processing/2024_11_18_whole_blood_replication_ali_et_all.tsv")

################################################################################
###       repeat in all tissues                                              ###
################################################################################

#+ run association test for all tissues using heteroplasmy genotype files ------
# init result object
res_list_het_calls <- list()
# loop through tissues --> we use the same tissue vector from above
# tissues <- unique(model_input_df$tissue)
# tissues <- tissues[-which(tissues == "kidney_cortex")]
for(t in seq_along(tissues)){
  # get tissue t
  tissue <- tissues[t]
  # tell us where we are
  cat(paste0("\n tissue is #",t,"/",length(tissues)))
  # load heteroplasmy genotypes
  t_genotype_file <- read_tsv(paste0("/genotype_files/split_by_tissue/",tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)
  t_genotype_file <- t_genotype_file %>% filter(Pos %in% pos_gene_df$Pos)
  t_genotype_file <- t_genotype_file %>% filter(signal_type %in% c("non_heteroplasmy","heteroplasmy")) 
  # drop feature_id in heteroplasmy file to avoid confusion
  t_genotype_file <- t_genotype_file %>% select(-feature_id)
  # add in tissue gene expression file
  t_gene_expression <- read_tsv(paste0("/gene_expression/",tissue,"_mtgenes_poly_A_peer_factors_regressed_out.txt"), show_col_types = FALSE)
  t_gene_expression <- melt(setDT(t_gene_expression), id.vars = "feature_id", variable.name = "SUBJID", value.name = "logTPM")
  # join tissue gene expression file to model input data
  t_df <-left_join(t_genotype_file,t_gene_expression, by = "SUBJID",relationship = "many-to-many") 
  # only consider the 7 position feature pairs where poly_A mRNA is immediately
  # upstream of p9 methylation site
  t_df <- t_df %>% filter(paste(Pos, feature_id) %in% paste(pos_gene_df$Pos, pos_gene_df$feature_id))
  # loop through positions
  positions <- unique(t_df$Pos)
  p_res_list_het_calls <- list()
  for(p in seq_along(positions)){
    # get position p
    position <- positions[p]
    # subset to position
    p_df <- t_df[which(t_df$Pos == position), ]
    # loop through features
    features <- unique(p_df$feature_id)
    f_res_list_het_calls <- list()
    for(f in seq_along(features)){
      # get feature f
      feature <- features[f]
      # subset to feature
      f_df <- p_df[which(p_df$feature_id == feature), ]
      # run association test
      single_sex_tissues <- c("ovary","prostate","testis","uterus","vagina")
      if(unique(p_df$tissue) %in% single_sex_tissues){   
        lm_m0 <- lm(sum_heteroplasmic_level ~ AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
        lm_m1 <- lm(sum_heteroplasmic_level ~ logTPM + AGE + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
      } else {
        lm_m0 <- lm(sum_heteroplasmic_level ~ AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
        lm_m1 <- lm(sum_heteroplasmic_level ~ logTPM + AGE + SEX + PC_1 + PC_2 + PC_3 + PC_4 + PC_5, data = f_df)
      }
      # get results and do statistical test
      res_m1 <- summary(lm_m1)
      LR_test <- lrtest(lm_m1,lm_m0) 
      # build res table
      tmp_df <- data.frame(tissue = tissue,
                           Pos = position,
                           feature_id = feature,
                           n_samples_het_calls = length(unique(p_df$SUBJID)),
                           property_tested = "heteroplasmy_calls",
                           estimate_het_calls =  lm_m1$coefficients["logTPM"],
                           lm_std_error_het_calls = res_m1$coefficients["logTPM","Std. Error"],
                           LR_pval_het_calls =  LR_test[2,  'Pr(>Chisq)'])
      
      # return result  
      f_res_list_het_calls[[p]] <- tmp_df
    }
    # return result 
    f_res_het_calls_df <- bind_rows(f_res_list_het_calls)
    p_res_list_het_calls[[p]] <- f_res_het_calls_df
  }
  # return result
  p_res_het_calls_df <- bind_rows(p_res_list_het_calls)
  res_list_het_calls[[t]] <- p_res_het_calls_df
}
# collapse list to data frame
res_het_calls_df <- bind_rows(res_list_het_calls)
# multiple testing correction: 
## study wide
res_het_calls_df$bonferroni_het_calls <- p.adjust(res_het_calls_df$LR_pval_het_calls,method = "bonferroni") 
res_het_calls_df$BH_het_calls <- p.adjust(res_het_calls_df$LR_pval_het_calls,method = "BH") 
## within tissue (as this seems closest to what Ali et al. had done.)
res_het_calls_df <- res_het_calls_df %>%
  group_by(tissue) %>%
  mutate(t_bonferroni_het_calls = p.adjust(LR_pval_het_calls, method = "bonferroni"),
         t_BH_het_calls = p.adjust(LR_pval_het_calls, method = "BH")) %>%
  ungroup()

## lets assign significance based on tissue wide bonferroni
res_het_calls_df <- res_het_calls_df %>% 
  mutate(sig_het_calls = case_when(t_bonferroni_het_calls <= 0.05 ~ TRUE,TRUE ~ FALSE))


#+ compare our heteroplasmy level results with ali et al -----------------------
# compare to previous
compare_df <- 
  left_join(compare_df,res_het_calls_df %>% mutate(Pos = as.factor(Pos))) %>% 
  select(tissue,Pos,n_samples_read_counts,n_samples_het_calls,gene_name,feature_id,poly_A,
         beta_gtex_ali_et_al,estimate_read_counts,estimate_het_calls,lm_std_error_read_counts,
         lm_std_error_het_calls,p_value_gtex_ali_et_al,LR_pval_read_counts,LR_pval_het_calls,
         t_bonferroni_read_counts,t_bonferroni_het_calls,sig_read_counts,sig_het_calls,sig_in_gtex_ali_et_al)

# join 2 result tables for all tissues and save the full table 
# for further analysis
both_tables_df <- 
  left_join(res_het_calls_df,
            res_df,
            by = c("tissue","Pos","feature_id"))
# save
write_tsv(both_tables_df,"~/data/mtDNA_variants/results/7_transcript_processing/2024_11_18_all_tissues_replication_ali_et_all.tsv")
