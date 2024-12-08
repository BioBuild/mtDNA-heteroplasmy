# author:  simon.wengert@helmholtz-muenchen.de
# purpose: test tissue specific enrichment for sequence alterations in various
#          genomic features. For example, formally test for CNS enrichment of
#          m1A/G methylations in tRNA transcript positions.
#

#+ source dependencies ---------------------------------------------------------
source("~/scripts/utils/global_settings.R")
library("ggbreak")
library("ggpubr")
mt_tRNA_modifications_df <- read_tsv("~/metadata/annotations/mt_tRNA_modified_sites_only.tsv")
m1A_G_methylations <- mt_tRNA_modifications_df %>% filter(rna_modification %in% c("m1A","m1G"))


#+ load annotated heteroplasmy table --------------------------------------------
gt_rna_var_df <- 
  readRDS(paste0(g_ops$metadata_dir,"2024_04_28_gt_rna_var_annotated.rds")) %>%
      mutate(molecular_process = case_when(Pos %in% m1A_G_methylations$genomic_pos ~ "rna_modification",
                                           TRUE ~ "dna_mutation")) %>%
  left_join(.,
            m1A_G_methylations %>%
              distinct(Pos = genomic_pos,nDNA_enzyme = confirmed_enzyme) %>%
              mutate(nDNA_enzyme = case_when(is.na(nDNA_enzyme) ~ "unknown", 
                                             TRUE ~ nDNA_enzyme)),
            by = "Pos")


#+ CNS enrichment testing rna_heteroplasmies -----------------------------------
# A CNS peak is clearly visible from the box-plot # heteroplasmies per tissue.
# The issue with CNS tissues in GTEx is that there are consistently fewer samples
# in CNS tissues and also that most of these samples have been obtained from the
# same donors. Therefore, (i) we are going to use a downsampling approach to formally
# test if the # of rna_heteroplasmies is in fact enriched in CNS tissues. On top
# of that (ii) we also have to take into account, that there are 3 nDNA encoded
# enzymes (one of them unknown). In fact, in the top 2.5 % wmh data there are 
# the follwoing # of methylation sites per nDNA_enzyme:
#
# # A tibble: 4 × 2
#  nDNA_enzyme     n
#  <chr>       <int>
#  TRMT10C        12
#  TRMT61B         1
#  unknown         2
#
# --> to make a fair comparison, we're also going to sub-sample the # of sites
# and pick only 1 site per enzyme. (actually, we're not going to do that for now)
#
# --> unit of comparison is the # of rna heteroplasmies per nDNA encoded RNA-
# methylation enzyme. 
#


#+ subset data to rna_heteroplasmies only for CNS enrichment testing -----------
test_for_enrichment_df <- 
  gt_rna_var_df %>%
    filter(molecular_process == "rna_modification") %>%  
    dplyr::select(SUBJID,tissue,tissue_category,Pos,sum_heteroplasmic_level,nDNA_enzyme)


#+ add N of donors -------------------------------------------------------------
test_for_enrichment_df <-
  left_join(test_for_enrichment_df,
            test_for_enrichment_df %>%
              distinct(SUBJID,tissue) %>%
              group_by(SUBJID) %>%
                summarise(N_tissues_per_donor = dplyr::n()) %>%
              ungroup(),
            by = 'SUBJID')


#+ apply CNS enrichment testing ------------------------------------------------
# set thr for random sample # tissues per donor and # sites per nDNA enzyme
n_tissues_to_sample <- 10
# n_sites_to_sample <- 1
# loop through iterations
n_iterations <- 1000
all_iterations <- 1:1000
res_list <- list()
for(n in seq_along(all_iterations)){
  
  cat(paste0("\niteration is #",n,"/",n_iterations))
  
  # 1st downsampling: # of donors per tissue
  #
  # select subset of individuals
  # require min n of tissues to be present for each donor and
  # sample the equal amount of tissues for that person
  N_tissue_subset_var_df <- 
    left_join(test_for_enrichment_df %>%
              filter(N_tissues_per_donor >= n_tissues_to_sample) %>% 
              group_by(SUBJID) %>%
                distinct(tissue) %>%
                # sample random tissues without replacement
                sample_n(size = n_tissues_to_sample,
                         replace = FALSE) %>%
              ungroup(),
              test_for_enrichment_df,
              by = c('SUBJID', 'tissue'))
  N_tissue_subset_var_df$is_cns <- ifelse(N_tissue_subset_var_df$tissue_category == 'CNS',TRUE,FALSE)
  
  # build contincency table & define nDNA enzyme as test feature 
  contingency_list <- list()
  test_features <- unique(N_tissue_subset_var_df$nDNA_enzyme)
  
  for(a in seq_along(test_features)){
    #  2nd downsampling: # of sites per nDNA enzyme
    # tmp_df <- N_tissue_subset_var_df %>%
    #             filter(nDNA_enzyme == test_features[a]) %>%
    #             group_by(tissue) %>%
    #              # sample random sites without replacement
    #               sample_n(size = n_sites_to_sample,
    #                        replace = FALSE) %>%
    #             ungroup()
  
    # build contingency table per feature to test
    tmp_df <- N_tissue_subset_var_df %>%
                    filter(nDNA_enzyme == test_features[a]) %>%
                    group_by(is_cns) %>%
                      summarise(!! paste0(test_features[a]) := dplyr::n()) %>%
                      #summarise(!! paste0(genomic_features[a]) := sum(N_modified_sites)) %>%
                    ungroup()
    tmp_df <- bind_cols(tmp_df,
                        N_tissue_subset_var_df %>%
                          filter(!nDNA_enzyme == test_features[a]) %>%
                          group_by(is_cns) %>%
                            summarise(!! paste0('non_',test_features[a]) := dplyr::n()) %>%
                            #summarise(!! paste0('non_',genomic_features[a]) := sum(N_modified_sites)) %>%
                          ungroup() %>%
                          dplyr::select(-is_cns))
    tmp_df <- as.data.frame(tmp_df %>% mutate(is_cns = c('non_cns','cns')))
    tmp_df <- tmp_df[c(2,1),]
    rownames(tmp_df) <- tmp_df[,1]
    tmp_df <- tmp_df[,-1]
    contingency_list[[a]] <- tmp_df
                    
  }
  names(contingency_list) <- test_features

  # apply Fisher's exact test 
  res_fisher_exact_list <- list()

  for(c in seq_along(contingency_list)){
      fisher_df <- fisher.test(contingency_list[[c]])
      res_fisher_exact_list[[c]] <- tibble(nDNA_enzyme = colnames(contingency_list[[c]])[1],
                                           odds_ratio = fisher_df$estimate,
                                           lower_conf_95 = fisher_df$conf.int[1],
                                           upper_conf_95 = fisher_df$conf.int[2],
                                           p_val = fisher_df$p.value) %>%
                                           mutate(ci_range = abs(upper_conf_95 - lower_conf_95))
    }
  res_fisher_exact_df <- do.call(rbind.data.frame, res_fisher_exact_list) 
  res_fisher_exact_df$iteration <- n 
  res_list[[n]] <- res_fisher_exact_df 
}
  
res_df <- do.call(rbind.data.frame, res_list) 

# save enrichment test results to disc for later use 
saveRDS(res_df, file = "2024_05_01_mtRNA_heteroplasmies_cns_enrichment_1000_iterations.rds")

#+ plot density curve n donors per tissue category -----------------------------
p_gt_density_n_donors_tissue_category <- 
    test_for_enrichment_df %>% 
        distinct(SUBJID,N_tissues_per_donor,tissue_category) %>%
        ggplot(aes(N_tissues_per_donor, colour = tissue_category)) +
        # geom_density(aes(y = ..scaled..)) +
        geom_density(size = 1.2) +
        geom_vline(xintercept = n_tissues_to_sample,
                   size = 1.5,
                   linetype = 'longdash') +
        scale_colour_manual(values = p_ops$colours_tissue_category) + #,
                            #guide = guide_legend(override.aes = list(colour = p_ops$colours_tissue_category, alpha = 1, shape = 15, size = 8))) +
        ylab('density') +
        xlab('# tissues per donor') +
        p_ops$my_theme +
        theme(legend.position = "bottom",
              legend.title = element_blank()) 

# save to disc
ggsave(p_gt_density_n_donors_tissue_category, 
       filename = "~/git/mtDNA_variants/paper/supplement/figures/2024_05_01_p_gt_density_n_donors_tissue_category.png",
       width =  8,
       height = 6,
       dpi = "retina")   

ggsave(p_gt_density_n_donors_tissue_category, 
       filename ="~/git/mtDNA_variants/paper/supplement/figures/2024_05_01_p_gt_density_n_donors_tissue_category.pdf",
       width =  8,
       height = 6,
       useDingbats = FALSE,
       dpi = 300)     


#+ calculate mean and confidence interval for enrichment tests -----------------

res_df <- 
  res_df %>% 
    group_by(nDNA_enzyme) %>%
      mutate(mean_OR = mean(odds_ratio),
             CI = 1.96*sd(odds_ratio)) %>%
    ungroup() %>%
    group_by(nDNA_enzyme) %>%
      mutate(lower_end = mean_OR - min(CI),
             upper_end = mean_OR + max(CI))


#+ add number of sites per nDNA enzymes for plotting ---------------------------
res_df <- 
  left_join(res_df,
            gt_rna_var_df %>%
              filter(molecular_process == "rna_modification") %>%
              distinct(Pos,nDNA_enzyme) %>%
              count(nDNA_enzyme) %>%
              mutate(n_methylation_sites = n),
            by = "nDNA_enzyme") %>%
    mutate(plot_label = paste0(nDNA_enzyme," (#",n_methylation_sites,")"))


#+ visualize enrichment test results -------------------------------------------
p_gt_cns_enrichment_n_of_heteroplasmies_1000_iterations <-
  res_df %>%
    ggplot(aes(mean_OR,plot_label, colour = nDNA_enzyme)) +
    geom_point(size = 3.5) +
    geom_vline(xintercept = 1) +
    xlab("odds ratio (#samples = 1000)") +
    geom_errorbar(aes(xmin = lower_end, xmax = upper_end),width=.2) +
    scale_x_continuous(breaks = round(seq(min(0.0), max(22.0), by = 1),1)) +
    scale_x_break(breaks = c(2,14), expand = TRUE, space = 0.3) +
    p_ops$my_theme +
    theme(axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.line.x.top = element_blank(),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12, face = "bold"),
          #legend.position = c(.82,.76),
           #legend.position = "none",
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 12)) +
    #guides(colour = guide_legend(title = "nDNA enzyme",
    #                             override.aes = list(shape = 15, size = 8))) + 
    # in case of all 3 enzymes in the data
    scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3"),
                       name = "nDNA enzyme" ,
                       guide = guide_legend(override.aes = list(shape = 15, size = 7, alpha = 1)))
    # in case of only 2
    #scale_color_manual(values = c("#d95f02","#7570b3"))

# save to disc
ggsave(p_gt_cns_enrichment_n_of_heteroplasmies_1000_iterations, 
       filename = "2024_05_01_p_gt_cns_enrichment_n_of_heteroplasmies_1000_iterations.png",
       #width =  5,
       #height = 3,
       width =  4,
       height = 3,
       dpi = "retina")   
ggsave(p_gt_cns_enrichment_n_of_heteroplasmies_1000_iterations, 
       filename ="2024_05_01_p_gt_cns_enrichment_n_of_heteroplasmies_1000_iterations.pdf",
       width =  5,
       height = 3,
       useDingbats = FALSE,
       dpi = 300)     


#+ generate polar coordinates plot annotation of m1A/G methylation sites -------
plot_methyl_sites_df <- 
  test_for_enrichment_df %>% 
    distinct(Pos,nDNA_enzyme) %>%
    mutate(fake_count = case_when(nDNA_enzyme == "TRMT10C" ~ 0.2,
                                  TRUE ~ 0.25))
# make plot
p_rna_heteroplasmy_sites_with_signal_all_tissues <-
  plot_methyl_sites_df %>%
    ggplot(aes(Pos,fake_count,label = Pos)) +
  #  geom_line(alpha = 0.6, size = 1.5) + 
    geom_point(aes(fill = nDNA_enzyme), size = 3.5, colour = "white",pch = 21) +
               #position = position_jitterdodge(jitter.height = 0.025)) +
    coord_polar(direction = -1) +
    scale_fill_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
    #scale_fill_manual(values = c("#d95f02","#7570b3")) +
    xlim(0,16569) +
    ylim(0,1) +
    theme_void() +
    theme(legend.position = "none") +
    geom_text_repel(nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.angle = 20) +
    annotate("text", x = 0, y = 0, label = "mtDNA\nposition", 
             size = 3.5, fontface = "bold", color = "black")  

# save to disc
ggsave(p_rna_heteroplasmy_sites_with_signal_all_tissues, 
       filename = "2024_05_01_p_rna_heteroplasmy_sites_with_signal_all_tissues.png",
       width =  6,
       height = 6,
       dpi = "retina")  
 ggsave(p_rna_heteroplasmy_sites_with_signal_all_tissues, 
       filename ="2024_05_01_p_rna_heteroplasmy_sites_with_signal_all_tissues.pdf",
       width =  6,
       height = 6,
       useDingbats = FALSE,
       dpi = 300)     
 
