# author:     simon.wengert@helmholtz-munich.de
#+ source dependencies ---------------------------------------------------------
source("~/git/mtDNA_variants/scripts/utils/global_settings.R")
library(data.table)
gt_gPCs <- read_tsv("gt_v8_gPCs_europ.tsv")

#+ load data -------------------------------------------------------------------
gt_rna_var_df <- readRDS("2024_04_28_gt_rna_var_passing_het_filters_df.rds")

#+ filter for samples and tissues passing cohort filters -----------------------
gt_rna_var_all_filters_df <- readRDS("2024_04_28_gt_rna_var_annotated.rds"))
passing_all_filters_mask_df <- gt_rna_var_all_filters_df %>%
  dplyr::select(biospecimen_repository_sample_id,SUBJID,tissue,Pos)
gt_rna_var_df <- gt_rna_var_df %>% 
  filter(biospecimen_repository_sample_id %in% unique(passing_all_filters_mask_df$biospecimen_repository_sample_id),
         tissue %in% unique(passing_all_filters_mask_df$tissue))

#+ saving time: optional filter for most common and abundant heteroplasmies ----
# build_mode <- "common_ones"
build_mode <- "all_positions"
if(build_mode == "common_ones"){
  # tell us which set of positions we are going to use
  cat(paste0("The following set of heteorplamies is used to build the genotype files: ",build_mode))
  # apply the filter
  gt_rna_var_df <- 
    left_join(gt_rna_var_df,
              passing_all_filters_mask_df %>%
                distinct(biospecimen_repository_sample_id,tissue,Pos) %>%
                mutate(keep = TRUE),
              by = c("biospecimen_repository_sample_id","tissue","Pos")) %>%
    filter(keep == TRUE)
} else if (build_mode == "all_positions"){
  # # tell us which set of positions we are going to use
  cat(paste0("The following set of heteorplamies is used to build the genotype files: ",build_mode," (after strand_bias filtering)"))
} else {
  stop("please make explicit which kind of heteroplasmy positions you want to use for building the genoypte files.")
}


#+ build genotype matrices per tissue ------------------------------------------
if(!file.exists(tmp_filename)){
  # init tissues for looping
  tissues <- unique(gt_rna_var_df$tissue)
  # loop for building long table genotypes
  t_list <- list()
  # track time
  start_time <- Sys.time()
  for(i in seq_along(tissues)){ 
      # select tissues
      tissue <- tissues[i]
      cat(paste0("tissue is: ",tissue," (# ",i,"/",length(tissues),")\n"))
      
      # subset for tissue in i
  	  gt_mutserve_tissue_df <- 
          gt_rna_var_df %>% 
              filter(tissue == tissues[i]) %>%
              mutate(feature_id = paste0('mt.',Pos))
      # now build genotype matrices only for the mtDNA positions passing wmh filter
      filtered_mtDNA_pos_ids <- gt_mutserve_tissue_df %>% distinct(feature_id) %>% pull(feature_id) 
      # this refers to the positions with confirmed coverage but no heteroplasmies
      gt_signal_detected_tissue_df <- 
          gt_rna_raw_df %>%
              filter(tissue == tissues[i],
                     Pos %in% gt_mutserve_tissue_df$Pos)
      gt_ldf <- list()
  	  gt_donors <- unique(gt_mutserve_tissue_df$SUBJID)
      # build genotype info from the ground up for each donor
      for(a in seq_along(gt_donors)){
              # build heteroplasmy-genotype matrix for
              # donor X from the ground up
              donor <- gt_donors[a]
              # add feature ids based on mtDNA positions which have signal - PER TISSUE!
              signal_detected <- 
                  gt_signal_detected_tissue_df %>%
                      filter(SUBJID == donor) %>%
                      mutate(feature_id = paste0('mt.',Pos))
              all_positions_df <- tibble(feature_id = sprintf("mt.%s",seq(1:16569)))
              # new - restrict to weighted median filtered positions
              all_positions_df <- 
                  all_positions_df %>% 
                      filter(feature_id %in% filtered_mtDNA_pos_ids)
              # get ids of mtDNA server variants
              variant_ids <- gt_mutserve_tissue_df %>%
                                  filter(SUBJID == donor) %>%
                                  distinct(feature_id) %>%
                                  pull(feature_id)
              # get ids of covered mtDNA positions - but not variants
              zero_ids <- signal_detected %>%
                              filter(!feature_id %in% variant_ids) %>%
                              distinct(feature_id) %>%
                              pull(feature_id)
              # get ids of mtDNA positions which don't have signal
              na_ids <- all_positions_df %>%
                              filter(!feature_id %in% c(variant_ids,zero_ids)) %>%
                              distinct(feature_id) %>%
                              pull(feature_id)
              # sanity check if the numbers match up
              stopifnot((length(na_ids)+length(zero_ids)+length(variant_ids)) == length(all_positions_df$feature_id)) 
              # assemble data frames
              variant_signal_df <- 
                  gt_mutserve_tissue_df %>%
                      filter(SUBJID == donor,
                             feature_id %in% variant_ids) %>% 
                      distinct(SUBJID,sum_heteroplasmic_level,Coverage,feature_id) %>%
                      filter(SUBJID == donor) %>%
                      # change data type for compatibility with other frames
                      mutate(sum_heteroplasmic_level =  as.character(round(sum_heteroplasmic_level, digits = 3)),
                            #sum_heteroplasmic_level = as.character(1), 
                            Coverage = as.character(Coverage),
                            signal_type = "heteroplasmy")
              zero_signal_df <- 
                  signal_detected %>%
                      filter(feature_id %in% zero_ids) %>%
                      # use real values for sum_heteroplasmic_level for positions with
                      # no variants called but sufficient coverage for being recorded
                      mutate(sum_heteroplasmic_level = as.character(round(sum_heteroplasmic_level, digits = 3)),
                             #sum_heteroplasmic_level = as.character(0), 
                             Coverage = as.character(sum_coverage),
                             signal_type = "non_heteroplasmy") %>%
                      distinct(SUBJID,feature_id,sum_heteroplasmic_level,Coverage,signal_type)		  
              no_signal_df <- 
                  tibble(SUBJID = donor,
                      feature_id = na_ids,
                      # use limix NaN encoding for missing valus
                      sum_heteroplasmic_level = "NaN",
                      Coverage = "NaN",
                      signal_type = "missing_value")
              # combine all data frames
              gt_genotypes_tmp_df <-
                  bind_rows(variant_signal_df,
                            zero_signal_df,
                            no_signal_df)
              # sort feature_id ascendingly in corrent order
              idx <- match(str_sort(gt_genotypes_tmp_df$feature_id, numeric = TRUE),gt_genotypes_tmp_df$feature_id)
              gt_genotypes_tmp_df <- gt_genotypes_tmp_df[idx, ]
              # store in list for further processing
              gt_ldf[[a]] <- gt_genotypes_tmp_df %>% distinct()
              # close donor loop
      }
      # collapse list to data frame
  	  gt_tissue_genotypes_df <- rbindlist(gt_ldf)
  	  # adding AD and BP columns for BB regression modelling 
      # (N_heteroplasmy counts vs. remaining counts)
      gt_tissue_genotypes_df <- 
    	  gt_tissue_genotypes_df %>%      
          mutate(N_modified_sites = round(as.double(sum_heteroplasmic_level) * as.integer(Coverage)),
                 AD = N_modified_sites,
                 BP = as.integer(Coverage) - N_modified_sites)
      # join with sample ids
      gt_tissue_genotypes_df <- 
          left_join(gt_tissue_genotypes_df, 
                    gt_rna_var_df %>% 
                      # use !! to make fitlering with the value
                      # stored in a variable work for dplyr::filter()
                      filter(tissue == !!tissue) %>%
                      distinct(SUBJID,biospecimen_repository_sample_id),
                    by = c('SUBJID'))
      
      # keep track of tissue colum
      gt_tissue_genotypes_df$tissue <- tissue
      # store for later
      t_list[[i]] <- gt_tissue_genotypes_df 
  }
  # print runtime
  end_time <- Sys.time()
  print(paste0("Run time ",round(abs(difftime(start_time, end_time, units = "mins",)),3)," min."))
  # [1] "Run time 7.501 min."
  # [1] "Run time 7.391 min."
  # collapse list to data frame
  gt_genotypes_df <- bind_rows(t_list)
  #
  # save a copy pre- genotype filtering
  saveRDS(gt_genotypes_df,tmp_filename)
  }else{
  # load existing file
  gt_genotypes_df <- readRDS(tmp_filename)
}

#+ post build filter for all pos or most frequent and variable het pos ---------
# selection_mode <- "frequent_and_variable_het_pos"
selection_mode <- "common_pos"
# selection_mode <- "all_positions"
if(selection_mode == "common_pos"){
  # tell us which set of positions we are going to use
  cat(paste0("Genotype files are subsetted to : ",selection_mode))
  # and now let's to posthoc filtering of genotypes built for all positions
  ## get Pos column back for merging
  tissue_pos_filter_mask_df <- passing_all_filters_mask_df %>% 
    mutate(feature_id = paste0("mt.",Pos)) %>%
    distinct(tissue,feature_id)
  ## merge and fitler
  # dimensions before: [1] 12379161       10
  gt_genotypes_df <- 
    left_join(gt_genotypes_df,
              tissue_pos_filter_mask_df %>%
                mutate(keep = TRUE),
              by = c("tissue","feature_id")) %>%
    filter(keep == TRUE) %>%
    dplyr::select(-c("keep"))
  # dimensions after: [1] 260196     10
} else if (selection_mode == "all_positions"){
  # # tell us which set of positions we are going to use
  cat(paste0("Genotype files are not filteres and the follwing het pos are used: ",selection_mode," (after strand_bias filtering)"))
} else {
  stop("please make explicit which kind of heteroplasmy positions you want to use keep in your genoypte files.")
}

#+ add cohort metadata and PCs -------------------------------------------------     
# join with genomic PCs
gt_genotypes_df <- left_join(gt_genotypes_df,gt_gPCs, by = c('SUBJID'))
# add cohort metadata    
gt_genotypes_df <- 
    left_join(gt_genotypes_df %>%
                mutate(tmp = feature_id) %>%
                separate(tmp,c("tmp_1","Pos"),sep = "mt.") %>%
                    dplyr::select(-tmp_1),
              gt_v8_lookup_df %>%
                    dplyr::select(-tissue),
              by = c('SUBJID', 'biospecimen_repository_sample_id')) %>%
    dplyr::select(feature_id,SUBJID,Pos,everything())

# add n_samples and n_positions per tissue
gt_genotypes_df <- 
  left_join(gt_genotypes_df,
            gt_genotypes_df %>% 
              distinct(tissue,biospecimen_repository_sample_id) %>%
              group_by(tissue) %>%
                summarise(n_samples = dplyr::n()) %>%
              ungroup()) %>%
  left_join(.,
            gt_genotypes_df %>% 
              distinct(tissue,Pos) %>%
              group_by(tissue) %>%
                summarise(n_positions = dplyr::n()) %>%
              ungroup())


#+ apply genotype based cohort_filters -----------------------------------------
# add columns for fraction of missing values per donor 
frac_signal_types_per_donor_df <- 
  gt_genotypes_df %>%
    group_by(SUBJID) %>%
      mutate(n_total = dplyr::n()) %>%
    ungroup() %>%
    group_by(SUBJID,signal_type,n_total) %>%
    count(SUBJID) %>%
    ungroup() %>%
    mutate(frac_type_per_donor = n/n_total) %>%
    dplyr::select(-c("n_total","n"))
# and per sample
frac_signal_types_per_sample_df <- 
  gt_genotypes_df %>%
    group_by(biospecimen_repository_sample_id) %>%
      mutate(n_total = dplyr::n()) %>%
    ungroup() %>%
    group_by(biospecimen_repository_sample_id,signal_type,n_total) %>%
    count(biospecimen_repository_sample_id) %>%
    ungroup() %>%
    mutate(frac_type_per_sample = n/n_total) %>%
    dplyr::select(-c("n_total","n"))
# join them to genotypes for filtering
gt_genotypes_df <- 
  left_join(gt_genotypes_df, 
            frac_signal_types_per_donor_df, 
            by = c("SUBJID", "signal_type")) %>%
  left_join(.,
            frac_signal_types_per_sample_df,
            by = c("biospecimen_repository_sample_id", "signal_type"))

#+ apply filters ---------------------------------------------------------------
gt_genotypes_df <- 
  gt_genotypes_df %>%
    # filter for donor missingness
    mutate(donors_to_filter_only = case_when(signal_type == "missing_value" ~ frac_type_per_donor,
                                             TRUE ~ 0.00)) %>%
    filter(donors_to_filter_only <= f_ops$donor_fraction_missing_sites) %>%
    # filter for sample missingness
    mutate(samples_to_filter_only = case_when(signal_type == "missing_value" ~ frac_type_per_sample,
                                              TRUE ~ 0.00)) %>%
    filter(samples_to_filter_only <= f_ops$donor_fraction_missing_sites) %>%
    # get rid of aux columns
    dplyr::select(-c(donors_to_filter_only,samples_to_filter_only))

#+ add postion variance column for filtering at testing stage ------------------
gt_genotypes_df <- 
  gt_genotypes_df %>%
    mutate(sum_heteroplasmic_level = as.double(sum_heteroplasmic_level)) %>%         
    group_by(feature_id) %>%
      mutate(var_per_pos = var(sum_heteroplasmic_level, na.rm = TRUE)) %>%
    ungroup()

#+ save results to disc --------------------------------------------------------
# one single file for all tissues
saveRDS(gt_genotypes_df, paste0("2024_04_28_gt_genotypes_all_t_",selection_mode,"_post_filtering.rds"))
# same content split by tissue for doing association tests per tissue
tissues <- unique(gt_genotypes_df$tissue)
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("tissue is: ",tissue," (#",t,"/",length(tissues),")\n"))
  t_df <- gt_genotypes_df %>% filter(tissue == !!tissue)
  # create directoy if doesn't exist
  save_dir <- paste0(g_ops$results_dir,"/split_by_tissue/",tissue,"/")
	ifelse(!dir.exists(file.path(save_dir)), 
			dir.create(file.path(save_dir), recursive = TRUE), 
			FALSE)
  # save to disc
  write_tsv(t_df, file = paste0('2024_04_28_gt_',tissue,"_",selection_mode,
           '_heteroplasmy_genotypes_long_format.tsv'))
}    



