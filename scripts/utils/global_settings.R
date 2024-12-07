# author: simon.wengert@helmholtz-muenchen.de
# purpose: define global settings and depencies for the project. This script 
#          should be sourced if required using source('utils/global_settings.R')


#+ source dependencies ---------------------------------------------------------
packages <- 
  c('tidyverse',
    'ggpubr',
    'cowplot',
    'plyranges',
    'ggrepel',
    'data.table')
lapply(packages, require, character.only = TRUE)


#+ define global options -------------------------------------------------------
g_ops <- list()
g_ops$date <- paste0(gsub("-","_",Sys.Date()),"_")


#+ get filter ops for all filters from csv file --------------------------------
all_filters_df <- read_csv("~/git/mtDNA_variants/metadata/utils/all_filters.csv")
# make filter options list for direct access in pipeline
f_ops <- as.list(all_filters_df$value)
names(f_ops) <- all_filters_df$name
f_ops$sites_seq_artefacts <- as.integer(unlist(strsplit(f_ops$sites_seq_artefacts, " ")))
f_ops$site_signal_type_exploration <- unlist(strsplit(f_ops$site_signal_type_exploration, split = " "))
f_ops$variant_type <- as.integer(f_ops$variant_type)
# 2024-06-24: we don't do this any more
# f_ops$top_wmh <- as.double(f_ops$top_wmh)
f_ops$min_cov_fwd <- as.integer(f_ops$min_cov_fwd)
f_ops$min_cov_rev <- as.integer(f_ops$min_cov_rev)
f_ops$delta_frac <- as.double(f_ops$delta_frac)
f_ops$fraction_sum_of_levels_per_pos <- as.double(f_ops$fraction_sum_of_levels_per_pos)
f_ops$site_variance_larger_than <- as.integer(f_ops$site_variance_larger_than)
f_ops$het_n_donors_min <- as.integer(f_ops$het_n_donors_min)
# 2024-06-24: we don't do this any more
# f_ops$mad_scaled_median <- as.double(f_ops$mad_scaled_median)
f_ops$tissue_n_donors_min <- as.integer(f_ops$tissue_n_donors_min)
f_ops$donor_fraction_missing_sites <- as.double(f_ops$donor_fraction_missing_sites)
f_ops$sample_fraction_missing_sites <- as.double(f_ops$sample_fraction_missing_sites)
f_ops$site_n_donors_min <- as.integer(f_ops$site_n_donors_min)


#+ define plotting option ------------------------------------------------------
p_ops <- list() 
p_ops$my_theme <- theme_classic() +
                  theme(plot.title = element_text(size=22, face = "bold"),
                        legend.title = element_text(size = 12, face = "bold"),
                        legend.text = element_text(size = 10, face = "bold"),
                        axis.text.x = element_text(size = 10, face = "bold"),
                        axis.text.y = element_text(size = 10, face = "bold"),  
                        axis.title.x = element_text(size = 14, face = "bold"),
                        axis.title.y = element_text(size = 14, angle = 90,face = "bold"),
                        strip.text = element_text(size = 12, face = "bold")) 
p_ops$grid_theme <- theme_bw() +
                  theme(plot.title = element_text(size=22, face = "bold"),
                        legend.title = element_text(size = 12, face = "bold"),
                        legend.text = element_text(size = 10, face = "bold"),
                        axis.text.x = element_text(size = 10, face = "bold"),
                        axis.text.y = element_text(size = 10, face = "bold"),  
                        axis.title.x = element_text(size = 14, face = "bold"),
                        axis.title.y = element_text(size = 14, angle = 90,face = "bold"),
                        strip.text = element_text(size=12, face = "bold")) 
p_ops$heatmap_theme <- 
  p_ops$my_theme + 
  theme(plot.title = element_text(size = 20, face = "bold"),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 22, face = "bold"), 
        legend.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1.5,"cm"))

# init custom colour scheme for tissue_categories
gt_v8_lookup_df <- read_tsv("~/git/mtDNA_variants/metadata/annotations/gt_v8_lookup_df.tsv")
tissue_categories <- unique(gt_v8_lookup_df$tissue_category)
colours_tissue_category <- c( "#FB8072", "#E1BE6A", "#D9D9D9",  "#B3DE69","#E66100", "#BEBADA",  "#8DD3C7", "#80B1D3",  "#FCCDE5" )
names(colours_tissue_category) <- tissue_categories
p_ops$colours_tissue_category <- colours_tissue_category

# init custom colour scheme for molecular event
p_ops$colours_molecular_event <-  c("#0072B2","#D55E00")
names(p_ops$colours_molecular_event) <- c("dna_mutation","rna_modification")

# create custom colour scale for mitochondrial complexes
# or let's just go with color brewer dark 2
# library(RColorBrewer)
# brewer.pal(n = 8,"Dark2")
p_ops$mt_complex_colours <- c("#1B9E77","#D95F02","#E6AB02","#E7298A", "#66A61E",  "#7570B3", "#666666")
names(p_ops$mt_complex_colours) <- c( "rRNA", "tRNA","complex_I", "complex_III", "complex_IV", "complex_V", "non_coding")
# create a custom color scale for codality (== MT transcript_type)
# old pink: "#CC79A7"
p_ops$colours_transcript_type <-
  c("#000000","#439894","#bf812d","#addd8e")
names(p_ops$colours_transcript_type) <- c("non_coding","protein_coding","rRNA_coding","tRNA_coding")

p_ops$colours_transcript_type_verbose <-
  c("#000000","#CC79A7","#bf812d","#addd8e", "#8c510a",
    "#c51b7d", "#01665e", "#35978f", "#80cdc1", "#c7eae5")
names(p_ops$colours_transcript_type_verbose) <- 
  c("non_coding","protein_coding","rRNA_coding","tRNA_coding","rRNA_coding_pos_947",
    "protein_coding_pos_834", "tRNA_coding_pos_9", "tRNA_coding_pos_16",
    "tRNA_coding_pos_37", "tRNA_coding_pos_58")

# adding positions: annotations rectangles (code from Na)
regions=read.table("~/git/mtDNA_variants/metadata/annotations/mtgenes.txt")
colnames(regions)=c("ID", "chr", "start", "end", "gene", "strand","mouse_gene","poly_A")
regions$colourf=as.character(regions$gene)
regions$colourf[grep("MT-T",regions$gene)]="MT-tRNA"
regions$colourf[grep("MT-RNR",regions$gene)]="MT-rRNA"
regions$colourf[grep("D-LOOP",regions$gene)]="D-LOOP"
library(RColorBrewer)
n <- length(unique(regions$colourf))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[8:(n+7)]
colref=data.frame(unique(regions$colourf),col_vector)
colnames(colref)=c("colourf", "colour")
colref$colourf=as.character(colref$colourf)
colref$colour=as.character(colref$colour)
regions$colour=colref$colour[match(regions$colourf,colref$colourf)]
regions = regions[order(regions$start),]
regions$colourf=as.factor(regions$colourf)
colref$colour=colref$colour[match(levels(regions$colourf),colref$colourf)]

protein_coding_genes <- 
  c("MT-ND1","MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
    "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6",
    "MT-CYB")
r_rna_coding_genes <- 
  c("MT-RNR1", "MT-RNR2")
t_rna_coding_genes <- 
  c("MT-TF", "MT-TV", "MT-TL1", "MT-TI", "MT-TQ", "MT-TM",
    "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-TS1",
    "MT-TD", "MT-TK", "MT-TG", "MT-TR", "MT-TH", "MT-TS2",
    "MT-TL2", "MT-TE", "MT-TT", "MT-TP")

regions$transcript_type <- 'non_coding'
regions$transcript_type <- ifelse(regions$gene %in% protein_coding_genes,
                                 'protein_coding',
                                  regions$transcript_type)
regions$transcript_type <- ifelse(regions$gene %in% r_rna_coding_genes,
                                  'rRNA_coding',
                                   regions$transcript_type)
regions$transcript_type <- ifelse(regions$gene %in% t_rna_coding_genes,
                                  'tRNA_coding',
                                  regions$transcript_type)
g_ops$mt_regions <- regions


#+ add annotation of single sex tissues ----------------------------------------
# g_ops$single_sex_tissues <- c"ovary","prostate","testis","uterus","vagina")


#+ defining global functions ---------------------------------------------------

# define function for extracting the trinucleotide profiles
get_trinucleotide_profiles <- function(data,fasta_path){

    # record run time
    start_time <- Sys.time()
    
    # dependencies
    library('Biostrings')

    # perform function
    fasta_file <- readDNAStringSet(fasta_path)
    data$ref_trinuc_profile <- "NA"

    positions <- data$Pos

    for(i in seq_along(positions)){
        
        minus_one <- positions[i]-1
        plus_one <- positions[i]+1
        data[i,"ref_trinuc_profile"] <- as.character(subseq(fasta_file, start=c(minus_one), end=c(plus_one)))
    
    }

    # calculate run time
    end_time <- Sys.time()
    print(paste0('run time: ',end_time - start_time, ' minutes'))

    return(data)
  }


#+ scratch ---------------------------------------------------------------------
if(F){
# usage
p <- 
  p + 
  geom_rect(data = regions, 
            inherit.aes = FALSE,
            colour = "black", 
            aes(xmin = start, 
                xmax = end, 
                ymin = -1, 
                ymax = -0.75, 
                fill = factor(regions$colourf)),
            alpha = 1)

p <- 
  p + 
  scale_fill_manual(values = colref$colour, 
                    name="mtDNA Loci")
}
