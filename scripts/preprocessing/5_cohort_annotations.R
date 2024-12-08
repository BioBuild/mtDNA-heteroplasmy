# author:   simon.wengert@helmholtz-munich.de
# notes:    add genomic and other annotations to cohort files 

#+ source dependencies ---------------------------------------------------------
source("~/scripts/utils/global_settings.R")

#+ load genuine heteroplasmy table ---------------------------------------------
gt_rna_var_df <- readRDS("2024_04_28_gt_rna_var_cohort_filters.rds")

  
#+ add info about energy consumption of tissues --------------------------------  
gt_rna_var_df <- 
  gt_rna_var_df %>%
      # add meta info about organ energy consumption  
      # see Wang et al 2010
      # https://academic.oup.com/ajcn/article/92/6/1369/4597507
      mutate(metabolic_activity = case_when(tissue_category == "cardiovascular" ~ 6,
                                            tissue_category == "CNS" ~ 5,
                                            tissue_category == "liver" ~ 4,
                                            tissue_category == "musculoskeletal_connective" ~ 3,
                                            tissue_category %in% c("others","digestive","blood_immune","endocrine","lung") ~ 2,
                                            tissue_category == "adipose" ~ 1)) 


#+ load annotations ------------------------------------------------------------
# RNA modified sites
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))
# genomic annotation
chrM_anno <- read_tsv("~/metadata/annotations/gencode_v35_annotation_chrM.gtf") %>% as_granges()


#+ annotate m1A/G rna methylations ---------------------------------------------
gt_rna_var_df <- 
  gt_rna_var_df %>%
       mutate(molecular_process = case_when(Pos %in% m1A_G_methylations$genomic_pos ~ "rna_modification",
                                            TRUE ~ "dna_mutation"))


#+ add genomic annotation ------------------------------------------------------     
stash_cols <- colnames(gt_rna_var_df)
gt_rna_var_df <-
  join_overlap_left(gt_rna_var_df %>% 
                       mutate(seqnames = 'chrM',
                              start = Pos,
                              end = Pos) %>%
                       as_granges(),
                     chrM_anno %>%
                       plyranges::select(gene_name)) %>%
                      as_tibble() %>%
                      dplyr::select(stash_cols,gene_name) %>%
                      distinct()


#+ save file to disc -----------------------------------------------------------
saveRDS(gt_rna_var_df,"2024_04_28_gt_rna_var_annotated.rds")

