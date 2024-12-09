# author:   simon.wengert@helmholtz-munich.de
# purpose:  load heteroplasmy genotype files, apply exact same filters as used 
#           in association testing, subset to RNA modifications and write out 
#           unique tissue-sample-pos lookup table for subsequent use in trans-
#           cript processing pipeline (based on pysam). The aim is to have the
#           read level summaries available for the same set of samples that went
#           into the association testing.

source("~/scripts/utils/global_settings.R")
data_path <- "/mtDNA_variants/"
model_input_path <- paste0(data_path,"metadata/genotype_files/split_by_tissue/")
tissues <- list.dirs(model_input_path, full.names = FALSE)[-1]
# load annotation of RNA modification sites 
mt_tRNA_modifications_df <- read_tsv("~/git/mtDNA_variants/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))

#+ load raw heteroplasmy estimates and allelic counts --------------------------
filename <- "2024_10_01_gt_rna_raw_position_lookup_table_rna_modifications.rds"
if(!file.exists(filename)){
  gt_rna_raw_df <- readRDS("2024_04_28_gt_raw_rna_heteroplasmy_estimates_allelic_counts_het_pos.rds")
  ## get sum coverage per Position
  gt_rna_raw_df <- gt_rna_raw_df %>%
    group_by(biospecimen_repository_sample_id,Pos) %>%
    mutate(sum_coverage = sum(Coverage)) %>%
    ungroup()
  ## subset to rna modifications
  gt_rna_raw_df <- gt_rna_raw_df %>% filter(Pos %in% m1A_G_methylations$genomic_pos)
  ## keep columns of interst
  gt_rna_raw_df <- gt_rna_raw_df %>% 
    distinct(Run,biospecimen_repository_sample_id,SUBJID,tissue,Pos,homoplasmic_base,sum_heteroplasmic_level,sum_coverage)
  ## save to disc
  saveRDS(gt_rna_raw_df,filename)
} else {
  gt_rna_raw_lookup_df <- readRDS(filename)
}

#+ loop through tissues and apply testing filters to each SNP ------------------
t_list <- list()
for(t in seq_along(tissues)){
  tissue <- tissues[t]
  cat(paste0("\n#tissues is ",t," out of #",length(tissues),"\n"))
  model_input_data <- read_tsv(paste0(model_input_path,tissue,"/2024_04_28_gt_",tissue,"_common_pos_heteroplasmy_genotypes_long_format.tsv"), show_col_types = FALSE)
  
  # only consider mtRNA modification sites 
  model_input_data <- model_input_data[which(model_input_data$Pos %in% m1A_G_methylations$genomic_pos), ] 
  
  ################################################################################
  ###                           filters                                        ###
  ################################################################################

  #+ drop NAs for now doing complete case analysis -------------------------------
  model_input_data <- model_input_data %>% drop_na()

  #+ filter out positions 0 variance ---------------------------------------------
  model_input_data  <- 
    model_input_data %>%
    filter(var_per_pos > f_ops$site_variance_larger_than)
  
  #+ filter for minimum number of donors -----------------------------------------
  min_donors_pos <- 
    model_input_data  %>%
    group_by(Pos) %>%
    count(Pos) %>% 
    filter(n >= f_ops$site_n_donors_min) %>%
    ungroup() %>%  
    pull(Pos)
  
  # keep tissue table for selecting bam files later
  t_df <- model_input_data %>% 
    distinct(SUBJID,tissue,Pos,signal_type,sum_heteroplasmic_level,Coverage)
    
  # store in list for later use
  t_list[[t]] <- t_df
  
}
# collapse list to data frame
all_tissue_sample_pos_df <- bind_rows(t_list)
 
#+ add sample-tissue-pos level info including homoplasmic base -----------------
all_tissue_sample_pos_df <- 
  left_join(all_tissue_sample_pos_df,
            gt_rna_raw_lookup_df %>%
                distinct(Run,SUBJID,biospecimen_repository_sample_id,
                         tissue,Pos,homoplasmic_base),
            by = c("SUBJID","tissue","Pos"),relationship = "many-to-many") %>%
    distinct()

#+ save to disc ----------------------------------------------------------------
write_tsv(all_tissue_sample_pos_df,"2024_10_01_gt_sample_list_at_rna_modification_sites.tsv")

#+ summarise read categories per sample and generate tissue level outs ---------
# loop through tissues
t_res_list <- list()
tissues <- unique(all_sample_pos_df$tissue)
for(t in seq_along(tissues)){
    # subset to tissue
    tissue <- tissues[t]
    t_df <- fread(paste0(data_path,input_date,"_",tissue,"_read_level_info_all_pos.csv"))
    cat(paste0("\ntissue is: ",tissue," #",t,"/",length(tissues),"\n"))
    
    # looping through positions
    p_res_list <- list()
    positions <- unique(t_df$heteroplasmy_position_one_based)
    for(p in seq_along(positions)){
        # subset to one-based position
        position <- positions[p]
        cat(paste0("\nposition is: ",position," #",p,"/",length(positions),"\n"))
        p_df <- t_df[which(t_df$heteroplasmy_position_one_based == position), ]
        # getting mtDNA coding strand for modification site of interest 
        # for inferring directionality when calculating cut site
        feature_strand <- m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "strand_name" ]  
        gene_name <- m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "gene_name" ] 

        # looping through samples        
        s_res_list <- list()
        samples <- unique(p_df$biospecimen_repository_sample_id)
        for(s in seq_along(samples)){
            # subset to samples
            sample <- samples[s]
            s_df <- p_df[which(p_df$biospecimen_repository_sample_id == sample), ]
            
            # use tryCatch to handle potential errors due to file 
            # non existing or little reands passing.
            file_status <- "file_passing"  
            
            tryCatch({                      
            
                # filtering and recording some basic numbers
                n_reads_total <- length(unique(s_df$read_id))
                # drop reads for which we don't have mate information available
                s_df <- s_df %>% filter(!is.na(mate_start))
                n_mates_avail <- length(unique(s_df$read_id))
                # drop reads fro which template_length == 0 (likely mapping issues)
                # I've never actually seen this in our data. More of a precaution.
                s_df <- s_df %>% filter(!template_length == 0)
                n_tlen_non_zero <- length(unique(s_df$read_id)) 
                # how many tlen values are < 0, > 0
                n_tlen_negative <- nrow(s_df[s_df$template_length < 0, ])      
                n_tlen_positive <- nrow(s_df[s_df$template_length > 0, ])

                # define DNA template start & stop coordinates
                ## if tlen is negative, then the mate maps before the read
                s_df <- s_df %>% 
                    mutate(template_start = case_when(template_length < 0 ~ mate_start, 
                                                    TRUE ~ read_start),
                        template_end = case_when(template_length < 0 ~ read_end,
                                                    TRUE ~ read_start + template_length))
                
                # keep track of cDNA template length summaries
                median_template_start = median(s_df$template_start)
                median_template_end = median(s_df$template_end)
                median_abs_tlen = median(abs(s_df$template_length))                

                # counting reads if they have been modified and cut      
                # 1) INFER MODIFICATION STATUS: assigned modification status to each read
                s_df <- s_df %>% mutate(is_modified = case_when(aligned_nucleotide != homoplasmic_base ~ TRUE,TRUE ~ FALSE))            
            
                # 2) INFER IF A READ HAS BEEN CUT OR NOT (TAKING PYSAM 0-based into accoung)
                if(feature_strand == "H-Strand"){

                    # assuming cut site as feature start with respect to mtDNA strand
                    s_df$cut_site_five_prime <-  (m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "feature_start"] -1)
                    s_df$cut_site_three_prime <-  (m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "feature_end"] -1)
                    # specify if read  is cut or not (assumption, based on cDNA template mapping coordinates)                  

                    ## five prime cut
                    s_df$is_cut_five_prime <- ifelse(s_df$template_start >= s_df$cut_site_five_prime, TRUE, FALSE)
                    # three prime cut
                    s_df$is_cut_three_prime <- ifelse(s_df$template_end <= s_df$cut_site_five_prime, TRUE, FALSE)
                
                } else if(feature_strand == "L-Strand") {

                    # assuming cut site as feature end with respect to mtDNA strand
                    s_df$cut_site_five_prime <-  (m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "feature_end"] - 1)
                    s_df$cut_site_three_prime <-  (m1A_G_methylations[which(m1A_G_methylations$genomic_pos == position), "feature_start"] - 1)
                    # specify if read  is cut or not (assumption!)
                    ## five prime
                    s_df$is_cut_five_prime <- ifelse(s_df$template_end <= s_df$cut_site_five_prime, TRUE, FALSE)
                    ## three prime
                    s_df$is_cut_three_prime <- ifelse(s_df$template_start >= s_df$cut_site_three_prime, TRUE, FALSE)
                
                } else {
                    stop("plase provide mtDNA strand annotation for RNA transcript")
                }


                # generate sample level summaries of read categories
                # GET NUMBERS OF READS "PER SPECIES" (REFERRING TO THE POSSIBLE 
                # COMBINATIONS OF BEING CUT AND BEING METHYLATED)
                ## get numbers with respect to five prime cut site
                n_is_modified_T_cut_five_prime_T <- nrow(s_df[which(s_df$is_modified == TRUE & s_df$is_cut_five_prime == TRUE), ])
                n_is_modified_F_cut_five_prime_T <- nrow(s_df[which(s_df$is_modified == FALSE & s_df$is_cut_five_prime == TRUE), ])
                n_is_modified_T_cut_five_prime_F <- nrow(s_df[which(s_df$is_modified == TRUE & s_df$is_cut_five_prime == FALSE), ])
                n_is_modified_F_cut_five_prime_F <- nrow(s_df[which(s_df$is_modified == FALSE & s_df$is_cut_five_prime == FALSE), ])

                ## get numbers with respect to three prime cut site
                n_is_modified_T_cut_three_prime_T <- nrow(s_df[which(s_df$is_modified == TRUE & s_df$is_cut_three_prime == TRUE), ])
                n_is_modified_F_cut_three_prime_T <- nrow(s_df[which(s_df$is_modified == FALSE & s_df$is_cut_three_prime == TRUE), ])
                n_is_modified_T_cut_three_prime_F <- nrow(s_df[which(s_df$is_modified == TRUE & s_df$is_cut_three_prime == FALSE), ])
                n_is_modified_F_cut_three_prime_F <- nrow(s_df[which(s_df$is_modified == FALSE & s_df$is_cut_three_prime == FALSE), ])

                ## get numbers with respect to both cut sites
                n_is_modified_T_cut_five_prime_T_cut_three_prime_T <- nrow(s_df[s_df$is_modified == TRUE & s_df$is_cut_five_prime == TRUE & s_df$is_cut_three_prime == TRUE, ])
                n_is_modified_T_cut_five_prime_T_cut_three_prime_F <- nrow(s_df[s_df$is_modified == TRUE & s_df$is_cut_five_prime == TRUE & s_df$is_cut_three_prime == FALSE, ])
                n_is_modified_T_cut_five_prime_F_cut_three_prime_T <- nrow(s_df[s_df$is_modified == TRUE & s_df$is_cut_five_prime == FALSE & s_df$is_cut_three_prime == TRUE, ])
                n_is_modified_T_cut_five_prime_F_cut_three_prime_F <- nrow(s_df[s_df$is_modified == TRUE & s_df$is_cut_five_prime == FALSE & s_df$is_cut_three_prime == FALSE, ])
                n_is_modified_F_cut_five_prime_T_cut_three_prime_T <- nrow(s_df[s_df$is_modified == FALSE & s_df$is_cut_five_prime == TRUE & s_df$is_cut_three_prime == TRUE, ])
                n_is_modified_F_cut_five_prime_T_cut_three_prime_F <- nrow(s_df[s_df$is_modified == FALSE & s_df$is_cut_five_prime == TRUE & s_df$is_cut_three_prime == FALSE, ])
                n_is_modified_F_cut_five_prime_F_cut_three_prime_T <- nrow(s_df[s_df$is_modified == FALSE & s_df$is_cut_five_prime == FALSE & s_df$is_cut_three_prime == TRUE, ])
                n_is_modified_F_cut_five_prime_F_cut_three_prime_F <- nrow(s_df[s_df$is_modified == FALSE & s_df$is_cut_five_prime == FALSE & s_df$is_cut_three_prime == FALSE, ])


                # ASSEMBLE RESULTS
                # count read info summaries for reporting
                # (bwlos should be n_is_modified over n_total... so there is an error
                # in here currently. Since this takes very long to compute I simply
                # re-calculated it the levels in 2024_07_15_read_level_analysis.R)
                read_counts_heteroplasmic_level <- nrow(s_df[which(s_df$is_modified == TRUE), ])/nrow(s_df)
                fraction_a <- nrow(s_df[which(s_df$aligned_nucleotide == "A"), ])/nrow(s_df)
                fraction_t <- nrow(s_df[which(s_df$aligned_nucleotide == "T"), ])/nrow(s_df)
                fraction_g <- nrow(s_df[which(s_df$aligned_nucleotide == "G"), ])/nrow(s_df)
                fraction_c <- nrow(s_df[which(s_df$aligned_nucleotide == "C"), ])/nrow(s_df)
                # get number of reads modified and number of reads overlapping
                n_reads_modified <- nrow(s_df[which(s_df$is_modified == TRUE), ])

            }, error = function(e) {
            
                # In case of an error, assign the specific value to file_status
                file_status <- "file_empty_or_no_reads_passing"
            
            })

            # also account for the possibility that there is 
            # no reads passing the filters above
            if(nrow(s_df) == 0){
                file_status <- "file_empty_or_no_reads_passing"
            }

            # ASSEMBLE RESULT DATA FRAMES
            if(file_status == "file_passing"){

                # build result table containig all summaries at the sample level       
                s_res_df <-
                    data.frame(
                        # sample specs
                        biospecimen_repository_sample_id = sample,
                        filename = unique(s_df$filename),
                        file_status = "passing",
                        tissue = tissue,
                        Pos = position,
                        gene_name = gene_name,
                        feature_strand = feature_strand,
                        date_processed = as.character(unique(s_df$date_processed)),
                        n_samples = length(samples),
                        homoplasmic_base = unique(s_df$homoplasmic_base),
                        heteroplasmic_level_mtDNA_server = unique(s_df$heteroplasmic_level_mtDNA_server), 
                        read_counts_heteroplasmic_level = read_counts_heteroplasmic_level,
                        fraction_a = fraction_a,
                        fraction_t = fraction_t,
                        fraction_g = fraction_g,
                        fraction_c = fraction_c,
                        coverage_mtDNA_server = unique(s_df$coverage_mtDNA_server),                  
                        # read processing/filter numbers
                        n_reads_total = n_reads_total,
                        n_reads_mates_avail = n_mates_avail,
                        n_reads_tlen_non_zero = n_tlen_non_zero,
                        n_reads_tlen_negative =  n_tlen_negative,    
                        n_reads_tlen_positive =  n_tlen_positive,
                        n_reads_passing = length(unique(s_df$read_id)), 
                        n_reads_modified = n_reads_modified,
                        # template_length specs
                        median_template_start = median_template_start,
                        median_template_end = median_template_end,
                        median_abs_tlen = median_abs_tlen,
                        # numbers with respect to five prime cut site 
                        n_is_modified_T_cut_five_prime_T = n_is_modified_T_cut_five_prime_T,
                        n_is_modified_F_cut_five_prime_T = n_is_modified_F_cut_five_prime_T,
                        n_is_modified_T_cut_five_prime_F = n_is_modified_T_cut_five_prime_F,
                        n_is_modified_F_cut_five_prime_F = n_is_modified_F_cut_five_prime_F,
                        # numbers with respect to three prime cut site
                        n_is_modified_T_cut_three_prime_T = n_is_modified_T_cut_three_prime_T,
                        n_is_modified_F_cut_three_prime_T = n_is_modified_F_cut_three_prime_T,
                        n_is_modified_T_cut_three_prime_F =  n_is_modified_T_cut_three_prime_F,
                        n_is_modified_F_cut_three_prime_F = n_is_modified_F_cut_three_prime_F,
                        # numbers with respect to both cut sites
                        n_is_modified_T_cut_five_prime_T_cut_three_prime_T = n_is_modified_T_cut_five_prime_T_cut_three_prime_T, 
                        n_is_modified_T_cut_five_prime_T_cut_three_prime_F = n_is_modified_T_cut_five_prime_T_cut_three_prime_F,
                        n_is_modified_T_cut_five_prime_F_cut_three_prime_T = n_is_modified_T_cut_five_prime_F_cut_three_prime_T,
                        n_is_modified_T_cut_five_prime_F_cut_three_prime_F = n_is_modified_T_cut_five_prime_F_cut_three_prime_F,
                        n_is_modified_F_cut_five_prime_T_cut_three_prime_T = n_is_modified_F_cut_five_prime_T_cut_three_prime_T,
                        n_is_modified_F_cut_five_prime_T_cut_three_prime_F = n_is_modified_F_cut_five_prime_T_cut_three_prime_F,
                        n_is_modified_F_cut_five_prime_F_cut_three_prime_T = n_is_modified_F_cut_five_prime_F_cut_three_prime_T, 
                        n_is_modified_F_cut_five_prime_F_cut_three_prime_F = n_is_modified_F_cut_five_prime_F_cut_three_prime_F 
                        )
            }

            if(file_status == "file_empty_or_no_reads_passing"){
                # return NAs for samples with no reads        
                s_res_df <-
                    data.frame(
                        # sample specs
                        biospecimen_repository_sample_id = sample,
                        filename = NA,
                        tissue = tissue,
                        file_status = "non_passing",
                        Pos = position,
                        gene_name = gene_name,
                        feature_strand = feature_strand,
                        date_processed = NA,
                        n_samples = NA,
                        homoplasmic_base = NA,
                        heteroplasmic_level_mtDNA_server = NA, 
                        read_counts_heteroplasmic_level = NA,
                        fraction_a = NA,
                        fraction_t = NA,
                        fraction_g = NA,
                        fraction_c = NA,
                        coverage_mtDNA_server = NA,                  
                        # read processing/filter numbers
                        n_reads_total = NA,
                        n_reads_mates_avail = NA,
                        n_reads_tlen_non_zero = NA,
                        n_reads_tlen_negative = NA,    
                        n_reads_tlen_positive = NA,
                        n_reads_passing = NA, 
                        n_reads_modified = NA,
                        # template_length specs
                        median_template_start = NA,
                        median_template_end = NA,
                        median_abs_tlen = NA,
                        # numbers with respect to five prime cut site 
                        n_is_modified_T_cut_five_prime_T = NA,
                        n_is_modified_F_cut_five_prime_T = NA,
                        n_is_modified_T_cut_five_prime_F = NA,
                        n_is_modified_F_cut_five_prime_F = NA,
                        # numbers with respect to three prime cut site
                        n_is_modified_T_cut_three_prime_T = NA,
                        n_is_modified_F_cut_three_prime_T = NA,
                        n_is_modified_T_cut_three_prime_F = NA,
                        n_is_modified_F_cut_three_prime_F = NA,
                        # numbers with respect to both cut sites
                        n_is_modified_T_cut_five_prime_T_cut_three_prime_T = NA, 
                        n_is_modified_T_cut_five_prime_T_cut_three_prime_F = NA,
                        n_is_modified_T_cut_five_prime_F_cut_three_prime_T = NA,
                        n_is_modified_T_cut_five_prime_F_cut_three_prime_F = NA,
                        n_is_modified_F_cut_five_prime_T_cut_three_prime_T = NA,
                        n_is_modified_F_cut_five_prime_T_cut_three_prime_F = NA,
                        n_is_modified_F_cut_five_prime_F_cut_three_prime_T = NA,
                        n_is_modified_F_cut_five_prime_F_cut_three_prime_F = NA
                        )
            }
            
            # store result
            s_res_list[[s]] <- s_res_df
        }
        # collapse and store
        p_res_df <- bind_rows(s_res_list)
        p_res_list[[p]] <- p_res_df
    }
    # collaps and store
    t_res_df <- bind_rows(p_res_list)
    # save tissue data frame 
    fwrite(t_res_df,paste0("/transcript_processing/tissue_summaries/",output_date,"_",tissue,"_sample_summary_from_read_level_analysis.tsv"))
    # go on 
    t_res_list[[t]] <- t_res_df
    # tell us where we are
    cat(paste0("\ntissue ",tissue," complete ------------------------------\n"))
}
# collaps to final result
all_tissues_res_df <- bind_rows(t_res_list)


#+ save final result to disc ---------------------------------------------------
fwrite(all_tissues_res_df,
       paste0("/transcript_processing/tissue_summaries/",output_date,"_all_tissues_sample_summary_from_read_level_analysis.tsv"),
       sep = "\t")    # tell us where we are
cat("\nall code execution complete ---------------------------------\n")


