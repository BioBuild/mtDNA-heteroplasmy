# mtDNA-heteroplasmy

## Description

This repository contains the code required to reproduce the analysis presented in the manuscript **"Tissue-specific apparent mtDNA heteroplasmy and its relationship with ageing and mtDNA gene expression"** by Wengert et al., which is currently under revision. You can access our preprint [here](link-to-paper). 

Below, you can find a conceptual overview of the project and the analyses conducted for this study:

![](README_files/project_overview.png)


- more details/short summary of the approach and key findings.


## Code Execution & Usability

The scripts in this repository work as part of a larger analysis pipeline, but may require minor adjustments to run on different systems. File paths and environment settings are specific to the original machines, so you’ll need to modify these to match your own setup. While the code is flexible and designed to be adaptable, it is not intended to run "out of the box" and may need configuration changes (e.g., paths, dependencies). That said, the  pipeline is modular and reproducible with proper setup, and we encourage contributions to improve its general usability.


## Analysis pipeline

   * [Variant calling](#Variant-calling)

   * [Data processing](#Data-processing) 
      
   * [Descriptive analysis](#Descriptive-analysis)

   * [Model evaluation](#Model-evaluation)

   * [Donor age testing](#Donor-age-testing)

   * [mtDNA cis-eQTL](#mtDNA-cis-eQTL)
   
   * [Mediation analysis](#Mediation-analysis)

   * [Transcript processing analysis](#Transcript-processing-analysis)


### Variant calling

- brief explanation
- code:
    XX. scripts for data download.
    XX. scripts for mtDNA-server variant calling.

### Data processing

- brief explanation (i.e. this is referring to Figure 1 in the manuscript.)
- code:
    1. `2023_05_17_gt_remote_processing.R`
    2. `2024_04_23_merge_raw_rna_allelic_counts_with_mtDNA_server_outs.R`
    3. `2024_04_23_heteroplasmy_filters.R`
    4. `2023_07_06_cohort_filters.R`
    5. `2023_07_15_cohort_annotation.R`
    6. `2023_07_07_prep_heteroplasmy_genotype_files.R`
    7. `2024_06_20_prep_gtex_celltype_deconvolution_files_plus_qc.R`

### Descriptive analysis

- brief explanation
- code:
    1. `pipeline_plotting.R`
    2. `2023_07_15_cohorts_plotting.R`
    3. `enrichment_testing.R`

### Model evaluation

- brief explanation
- code: 
    1. `2024_04_28_derive_simulation_params.R`
    2. `2023_11_22_betabinomial_heteroplasmy_simulations.ipynb`
    3. `2023_12_13_model_assessment_simulations.R`

### Donor age testing

- brief explanation
- code:
    1. `2023_06_15_lm_analytical_all_tissues.R`
    2. `2024_24_28_LM_A_qc_plots.R`
    3. `2023_09_15_bb_permutations_per_tissue.R`
    4. `2023_07_14_BB_fit_analytilcal_pvals.R`
    5. `2023_07_14_BB_calc_empirical_pvals.R`
    6. `2023_08_07_BB_annotate_results.R`
    7. `2023_08_28_donor_age_model_diagnostics.R`
    8. `2024_06_19_donor_age_cell_type_interaction.R`
    9. `2023_07_15_plot_pheno_tests`

### mtDNA cis-eQTL

- brief explanation
- code: 
    1. `2024_06_13_regress_out_peer_factors_from_gene_expression.R`
    2. `2023_06_19_lm_cis_eQTL_analytical_all_tissues.R`
    3. `2024_24_28_cis_eQTL_LM_A_qc_plots.R`
    4. `2023_08_03_BB_permutations_per_tissue.R`
    5. `2023_08_08_BB_cis_eQTL_analytical_all_tissues.R`
    6. `2023_08_08_BB_calc_empirical_pvals.R`
    7. `2023_08_10_cis_eQTL_merge_and_annotate_results.R`
    8. `2023_08_10_cis_eQTL_model_diagnostics.R`
    10. `2024_06_17_ciseQTL_cell_type_interaction_analysis.R`
    11. `2024_06_18_cell_type_fractions_per_tissue.R`
    xx. `2024_08_14_sanity_check_file_ids.R`
    12. `2023_08_01_plot_mol_test.R`

### Mediation analysis

- brief explanation
- code:
    1. `2024_04_02_mediation_analysis.`
    2. `2024_09_23_dono_age_logTPM_interaction.R`

### Transcript processing analysis

- brief explanation
- code:
    1. `2024_10_01_select_tissue_sample_pos_bamfiles.R`
    2. `slurm_submit.py`
    3. `slurm_build_tissue_command.sh`
    4. `2024_09_25_count_reads_at_cutsites.py`
    5. `2024_09_28_sample_level_summary_reads_at_cutsites.R`
    6. `2024_08_23_qc_descriptive_sample_level_read_summaries.R`
    7. `2024_11_18_replication_analysis.R`
    8. `2024_11_18_position_level_analysis.R`
    9. `2024_11_18_sample_level_analysis.R`
    10. `2024_10_03_mt_transcript_lengths.R`


## Dependencies

- `python` packages:

- `R` packages:


## Contact

We are grateful for any feedback or questions about the analysis code! 

- **Coding-Related Questions:**  
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

