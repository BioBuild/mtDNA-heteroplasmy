# mtDNA-heteroplasmy

# Description

This repository contains the code required to reproduce the analysis presented in the manuscript **"Tissue-specific apparent mtDNA heteroplasmy and its relationship with ageing and mtDNA gene expression"** by Wengert et al., which is currently under revision. You can access our preprint [here](https://www.biorxiv.org/content/10.1101/2024.12.11.627989v1.abstract). 

Below, you can find a conceptual overview of the project and the analyses conducted for this study:

![](README_files/project_overview.png)


# Code usability and data access

The scripts in this repository work as part of a larger analysis pipeline, but may require minor adjustments to run on different systems. We also do not provide individual level data (either raw data or processed data as part of our project) in any part of our repository, these would have to be obtained through gaining [approved access with GTEx](https://gtexportal.org/home/protectedDataAccess). As such, scripts provided here are not intended to run as is but for demonstration of our analyses.

# Obtain mtDNA alignments from bulk RNAseq/WGS data and perform variant calling 

1. We obtain reads mapping to the rCRS mtDNA reference genome [NC_012920 ](https://www.ncbi.nlm.nih.gov/nuccore/251831106) (which is part of the human [hg38 reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)) using [samtools](https://www.htslib.org/), using flag [-F 3852](https://broadinstitute.github.io/picard/explain-flags.html) to select only the reads that do not fall under the following categories: read unmapped (0x4), mate unmapped (0x8), not primary alignment (0x100), read fails platform/vendor quality checks (0x200), read is PCR or optical duplicate (0x400), supplementary alignment (0x800). 
   
```
samtools view -b -F 3852 $file.bam -o $file.mtdna.bam chrM
samtools index $file.mtdna.bam
``` 
2. We perform variant calling using [mtDNA-server](https://github.com/seppinho/mutserve?tab=readme-ov-file) for each sample individually, then combine the variant calls to form our variant call set, using the [rCRS.fasta](https://raw.githubusercontent.com/seppinho/mutserve/master/files/rCRS.fasta) file as reference. 
```
reffile="rCRS.fasta"
java -jar mutserve-1.3.0.jar analyse-local --input $file.mtdna.bam --reference $reffile --level 0.01 --output $file.mtdna.var.txt 
```

# Variant data filtering and processing

We perform the following steps for apparent heteroplasmy variant filtering and processing, all scripts for these steps are in the ```scripts/preprocessing``` directory 

1. concatenating and preprocessing of raw variant call files from mtDNA-server, as shown in  ```1_variant_files_preprocessing.R```
2. identifying the likely inherited allele using the homoplasmic alleles in each individual at each apparent heteroplasmic position, using variant calls from WGS in Whole Blood, as shown in ```2_wgs_alleles.R```
3. applying quality control filters on apparent heteroplasmy calls, as shown in ```3_variant_filters.R```
4. apply cohort filters to identify common apparent heteroplasmies with adequate variance between individuals for association testing, as shown in ```4_cohort_filters.R```
5. add in cohort annotations, as shown in ```5_cohort_annotations.R```
6. prepare apparent heteorplasmy genotype files for association testing, as shown in ```6_prep_het_genotype_file.R```
7. identify cell types with high median [xCell scores](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) in [7 cell types in 35 GTEx tissues](https://www.science.org/doi/10.1126/science.aaz8528) for use in celltype interaction analyses, as shown in ```7_celltype_proportions.R```

# Landscape of apparent mtDNA heteroplasmy and enrichment testing in CNS tissues 

We perform the following explorations for the landscape of mtDNA heteroplasmy and mtRNA modifications in 49 tissues in GTEx v8, as shown in the ```scripts/landscape``` directory 

1. identification of numbers and median variant allele frequencies of mtDNA heteroplasmy and mtRNA modifications across donors in each tissue
2. identification of transitions and transversions in mtDNA heteroplasmy
3. identification of putative inherited/somatic nature mtDNA heteroplasmy through looking into the number of tissues they occur in as described in [An et al Nat Genet 2024](https://www.nature.com/articles/s41588-024-01838-z) 
4. obtaining associations between putative inherited/somatic nature mtDNA heteroplasmy and median mtDNA copy number (mtDNA-CN) as previously established in [Rath et al PNAS 2024](https://www.pnas.org/doi/10.1073/pnas.2402291121)
5. checking for tissue occurrence of pathogenic mtDNA heteroplasmy at mt.3243A>G as reviewed in [Gomes et al Hum Mol Genet 2021](https://doi.org/10.1093/hmg/ddab156), previously associated with various clinical phenotypes ranging from severe MELAS to mild deafness and glucose intolerance
6. enrichment testing of occurrence of mtRNA modifications in CNS tissues

1-5 are detailed in ```landscape_explorations.R``` and 6 is detailed in ```cns_enrichment.R```

# Evaluation of association testing models for RNAseq derived apparent heteroplasmy 

We obtain realistic relationships between donor age, RNAseq coverage on mtDNA genes, and apparent heteroplasmy using data from GTEx v8, and use these parameters to simulate a realistic null situation where apparent heteroplasmy is related to RNAseq coverage but not donor age. We then use this simulated data to test if a linear model (LM) or a beta-binomial model (BB) has well-calibrated false discovery rates (FDR). Scripts for this are in the ```scripts/modeltesting``` directory

1. obtain realistic parameters in ```get_realistic_params.R```
2. simulations in ```simulations.ipynb```
3. testing LM and BB models on the simulated data in ```test_simulations.R```

# Donor age association testing

From the evaluation of association testing models we determined that the BB model is better calibrated than the LM especially when RNAseq coverage is low and heteroskedasticity of heteorplasmy is high. However BB takes a much longer time than LM, and even more so with permutations that can allow us to derive empirical p values. As such in our paper we detail our recommendation for a two-step approach where we perform LM first then BB only on the significant associations. We use this two-step approach in the paper to examine the effect of donor age on apparent heteroplasmy VAF. Scripts for this are in the ```scripts/associations``` directory 

1. LM test for donor age associations with VAF in ```donorage_lmtest.R```
2. BB test functions for donor age associations with VAF in ```functions_bb_empirical.R```
3. BB test for donor age associtions with VAF in ```donorage_bbtest.R```
4. BB permutations and derivation of empirical p values are in ```donorage_bbpermute.R``` and ```donorage_bbempirical.R```

We also investigate the interaction effect of donor age and cell type [xCell scores](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) in [7 cell types in 35 GTEx tissues](https://www.science.org/doi/10.1126/science.aaz8528). This is shown in ```donorage_celltype_interactions.R```

# Apparent heteroplasmy eQTL  

We further use the same two-step approach in the paper to examine the effect of apparent heteroplasmy VAF on mtDNA gene expression, in polyA mtDNA genes (all protein coding genes except for _MT-ND6_). For all mtDNA polyA gene expressions we use residuals of their log(TPM+1) gene expression levels after correcting for tissue-specific PEER factors, as shown in ```get_geneexp_residuals.R```. Scripts for this are in the ```scripts/eqtl``` directory 

1. LM test for VAF eqtls in ```eqtl_lmtest.R```
3. BB test for VAF eqtls in ```eqtl_bbtest.R```
4. BB permutations and derivation of empirical p values are in ```eqtl_bbpermute.R``` and ```eqtl_bbempirical.R```

We also investigate the interaction effect of cell type [xCell scores](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) on VAF eqtls in [7 cell types in 35 GTEx tissues](https://www.science.org/doi/10.1126/science.aaz8528). This is shown in ```eqtl_celltype_interactions.R```

# Mediation analysis between apparent heteroplasmy, mtDNA gene expression and donor age 

We then perform two analyses to test the relationships between apparent heteroplasmy, mtDNA gene expression and donor age for 10 instances where the apparent heteroplasmy eqtl also has a donor age association. Scripts for both analyses are in ```scripts/mediation```

1. Mediation analysis using partial correlations in ```mediation_analysis.R```
2. interaction effects between donor age and mtDNA gene expression on apparent heteroplasmy VAF in ```donorage_eqtl_interaction.R```
   
# mt-tRNA modification roles in transcript processing 

Finally we perform a few analyses to test whether p9 mt-tRNA modifications affect the gene expression levels of the mtDNA genes on their 5' ends, as previously proposed in many studies including [Meynier et al Nat Comms 2024](https://www.nature.com/articles/s41467-024-49132-0), and tested implicitly in [Ali et al Comms Bio 2020](https://www.nature.com/articles/s42003-020-0879-3). Scripts for these analyses are in ```scripts/transcript```

1. Replication of results shown in [Ali et al Comms Bio 2020](https://www.nature.com/articles/s42003-020-0879-3) in Whole Blood, then performed in all other tissues in GTEx v8, in ```p9_5prime_replication.R```
2. Same analysis done only using modified mt-tRNA reads per sample, asking if a 5' cut on the same reads are associated with 5' gene expression, in ```p9_5prime_cut.R```
3. Read level data per sample across all tissues at p9 mt-tRNA modifications are obtained using [pysam](https://pysam.readthedocs.io/en/latest/api.html) as shown in ```p9_count_cutreads.py``` and annotated summarised in ```p9_read_level_data.R``` 

# Contact

We are grateful for any feedback or questions about the analysis code! 

- **Code-Related Questions:**  
  If you have questions or encounter issues with the code, please submit an issue via `github`.

- **Scientific Correspondence:**  
    For scientific correspondence please reach out directly to:

  - Simon Wengert: [simon.wengert@helmholtz-munich.de](mailto:simon.wengert@helmholtz-munich.de)  
  - Dr. Na Cai: [na.cai@helmholtz-munich.de](mailto:na.cai@helmholtz-munich.de)


## Licensing

This project is licensed under the MIT License. However, the authors respectfully request that it be used only for non-commercial purposes, unless prior written consent is obtained.

### What this means

- **Non-commercial use** includes academic research, educational projects, personal experimentation, and open collaboration.
- **Commercial use** includes, but is not limited to:
  - Use within a for-profit entity (e.g., private companies, corporate labs).
  - Redistribution of the software for monetary gain.
  - Embedding the software in a product or service that is sold or licensed.

If you wish to use this software for commercial purposes, please contact the authors to discuss licensing terms.

