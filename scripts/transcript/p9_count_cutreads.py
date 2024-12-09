# author:   simon.wengert@helmholtz-munich.de
# purpose:  pysam based implementation to extract all reads overlapping
#           rna modification sites and save these files to disc for heteroplasmy
#           paper transcript processing analysis. This script does the extraction
#           for one tissue.


#+ source dependencies ---------------------------------------------------------
import os
import sys
import pysam
import pandas as pd
import session_info


#+ params, specs and helper files ----------------------------------------------
# option parsing
tissue = sys.argv[1]    
date = sys.argv[2]      
save_path = sys.argv[3]
# define output name
output_filename = f"{save_path}/{date}_{tissue}_read_level_info_all_pos.csv"

# keep processing log
log_data = []
# account for differences in 0-based vs. 1-based genomic coordinates
offset = 1
# load tissue sample rna modfication combinations for selecting reads
all_sample_pos_df = pd.read_csv("2024_10_01_gt_sample_list_at_rna_modification_sites.tsv", sep="\t")
# data paths
gt_dir="/mtbams/"


# + define read and base level filter options ----------------------------------
filter_options = {
    'min_mapping_quality': 20,
    'min_base_quality': 20,
    'min_alignment_quality': 30,
    'exclude_secondary': True,
    'exclude_supplementary': True,
    'exclude_duplicate': True
}

#+ define function to extract read level information ---------------------------
def extract_read_level_info(filename, modified_position_zero_based, filter_options, biospecimen_repository_sample_id, date, tissue, pos, homoplasmic_base, heteroplasmic_level_mtDNA_server, coverage_mtDNA_server):
    # open the BAM file
    bamfile = pysam.AlignmentFile(filename, "rb")
    # create pysam's pileup iterator only at RNA modified position 
    try:
        pileup = bamfile.pileup("MT", modified_position_zero_based, modified_position_zero_based + 1)
    except ValueError:
        pileup = bamfile.pileup("chrM", modified_position_zero_based, modified_position_zero_based + 1)
    # init dictionary to contain reads in loop
    unique_reads = {}
    # initi set to keep track of processed read IDs for keeping unique reads only
    processed_reads = set()
    # iterate over pysam pileup columns
    for pileup_column in pileup:
        # only process reads at the modified genomic position (modified_position_zero_based)
        if pileup_column.pos == modified_position_zero_based:
            # iterate over individual reads
            for pileup_read in pileup_column.pileups:
                read = pileup_read.alignment
                alignment_read_query_pos = pileup_read.query_position
                read_id = read.query_name
                # apply read level filters
                if (
                    not read.is_unmapped  
                    and (not filter_options['exclude_secondary'] or not read.is_secondary)
                    and (not filter_options['exclude_supplementary'] or not read.is_supplementary)
                    and (not filter_options['exclude_duplicate'] or not read.is_duplicate)
                    and (read.mapping_quality >= filter_options['min_mapping_quality'])
                ):
                    # only process if the query position is valid
                    if alignment_read_query_pos is not None and isinstance(alignment_read_query_pos, int):
                        if 0 <= alignment_read_query_pos < len(read.query_sequence):
                            aligned_nucleotide = read.query_sequence[alignment_read_query_pos]
                            base_quality = read.query_qualities[alignment_read_query_pos]  
                            # check base quality and alignment quality per base
                            if base_quality >= filter_options['min_base_quality'] and read.mapping_quality >= filter_options['min_alignment_quality']:
                                # fill read_data result object
                                read_data = {
                                    'date_processed': date,
                                    'genomic_coordinates_default': "zero_based",
                                    'biospecimen_repository_sample_id': biospecimen_repository_sample_id,
                                    'filename': os.path.basename(filename),
                                    'tissue': tissue,
                                    'heteroplasmy_position_one_based': pos,
                                    'homoplasmic_base': homoplasmic_base,
                                    'aligned_nucleotide': aligned_nucleotide,
                                    'base_quality': base_quality,
                                    'heteroplasmic_level_mtDNA_server': heteroplasmic_level_mtDNA_server,
                                    'coverage_mtDNA_server': coverage_mtDNA_server,
                                    'modified_position_one_based': modified_position_one_based,
                                    'modified_position_zero_based': modified_position_zero_based,
                                    'read_id': read_id,
                                    'read_start': read.reference_start,
                                    'read_end': read.reference_end,
                                    'mapping_quality': read.mapping_quality,
                                    'cigar_string': read.cigarstring,
                                    'mate_start': read.next_reference_start if read.next_reference_start != -1 else "NA",
                                    'template_length': read.template_length 
                                }
                                # Append read level info to sample level dictionary
                                unique_reads[read_id] = read_data
    # collapse to Pandas DataFrame
    df = pd.DataFrame.from_records(list(unique_reads.values()))
    # return the DataFrame
    return df


#+ extract read level information per tissue -----------------------------------
# subset to tissue of interest
tissue_df = all_sample_pos_df[all_sample_pos_df['tissue'] == tissue]
# get the samples for the current tissue
samples = tissue_df['biospecimen_repository_sample_id'].unique()
# initialize a list to collect sample-level DataFrames for the current tissue
sample_dfs = []
# iterate over each sample_id in tissue_pos_df
for iterator, sample_id in enumerate(samples, 1):
    print(f"sample is: #{iterator}/{len(samples)}")  # current sample status
    filename = os.path.join(gt_dir, f"{sample_id}.mtfiltered.bam")  
    # subset to sample data frame      
    sample_df = tissue_df[tissue_df['biospecimen_repository_sample_id'] == sample_id]
    # loop through each position
    for pos in sample_df['Pos']:
        # specify sample level variables to transfer from original heteroplasmy data
        pos_row = sample_df[sample_df['Pos'] == pos].iloc[0] 
        homoplasmic_base = pos_row['homoplasmic_base']  
        heteroplasmic_level_mtDNA_server = pos_row['sum_heteroplasmic_level']  
        coverage_mtDNA_server = pos_row['Coverage']
        # clearly state zero/one based thing
        modified_position_one_based = pos # all our other analyses
        modified_position_zero_based = pos - offset # pysam
        # run extraction
        try:
            sample_pos_reads_df = extract_read_level_info(
                filename=filename, 
                modified_position_zero_based=modified_position_zero_based,  # Adjusting to zero-based
                filter_options=filter_options,  
                biospecimen_repository_sample_id=sample_id,  
                date=date, 
                tissue=tissue, 
                pos=pos, 
                homoplasmic_base=homoplasmic_base, 
                heteroplasmic_level_mtDNA_server=heteroplasmic_level_mtDNA_server, 
                coverage_mtDNA_server=coverage_mtDNA_server
            )
            command_successful = True  # Set success flag
            sample_dfs.append(sample_pos_reads_df)  # Append DataFrame for the current position
        except Exception as e:
            print(f"Error processing {sample_id} at position {pos}: {e}")
            command_successful = False
        # keep track of specs and call params
        log_data.append({
            'date': date,
            'sample_id': sample_id,
            'tissue': tissue,
            'pos': pos,
            'command_successful': command_successful
        })
    # collapse all sample-level DataFrames for the current tissue into one DataFrame
    if sample_dfs:  # Ensure there are DataFrames to concatenate
        tissue_combined_df = pd.concat(sample_dfs, ignore_index=True)
        # save to disc:
        ## result data frame
        tissue_combined_df.to_csv(output_filename, index=False)
        ## log file 
        ## save log_df to save_path as a tab delimited file
        log_df = pd.DataFrame(log_data)
        # save log_df to save_path as a tab delimited file
        log_df.to_csv(os.path.join(save_path, f'log_files/{date}_{tissue}_processing_log.csv'), sep="\t", index=False)


#+ session info ----------------------------------------------------------------
session_info.show()
