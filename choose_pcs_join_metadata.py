#!/usr/bin/env python3
# script to merge PCs from different input data types (e.g., genetic variants, regional methylation, gene expression, etc.) with metadata/covariates for downstream analysis

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Choose PCs from multiple input PC files and merge chosen PCs with metadata/covariates table.")
    # specify path to input metadata file
    parser.add_argument("--input_metadata", required=True, help="Path to initial input metadata/covariates file to be joined with PCs of interest.")
    # specify path to input pc files
    parser.add_argument("--input_pcs", required=True, nargs="+", help="Path to input PC files generated or preprocessed with make_pcs_stepwise.py.")
    # specify pc counts
    parser.add_argument("--pc_counts", required=True, nargs="+", help="Number of PCs starting with PC1 to retain for input PC files, in order of inputs specified for --input_pcs.")
    # specify prefix for output files
    parser.add_argument("--output_prefix", required=True, help="Specify prefix for output merged metadata/covariates/PCs file and covariate correlation heatmap.")
    # return arguments
    return parser.parse_args()
    
# join tables on sample columns
def join_metadata_pc_tables(input_metadata_data_frame,input_pc_data_frames,input_pc_counts):
    # set up merged data frame
    merged_metadata_pc_data_frame = input_metadata_data_frame
    # loop through input pc tables
    for idx, i in enumerate(input_pc_data_frames):
        # subset for first x pcs specified in corresponding element of input_pc_counts array
        current_input_pc_data_frame_subset=i.iloc[:,0:(input_pc_counts[idx]+1)]
        # merge this subset with initial data frame
        # merge subset pc table with existing merged data table on SAMPLE column using inner join
        merged_metadata_pc_data_frame=merged_metadata_pc_data_frame.merge(current_input_pc_data_frame_subset,on='SAMPLE')
    return(merged_metadata_pc_data_frame)
    
# make covariate/PC correlation heatmap to evaluate relationships between PCs and covariates included in downstream regressions
def make_covariate_correlation_heatmap(final_metadata_pc_data_frame,output_prefix):
    # separate covariates
    covariates = final_metadata_pc_data_frame.drop(columns=['SAMPLE'])
    # convert categorical variables into dummies
    non_numeric_cols = covariates.select_dtypes(include=['object', 'category']).columns
    dummy_vars = pd.get_dummies(covariates[non_numeric_cols], drop_first=True)
    # get numerical variables
    numeric_vars = covariates.drop(columns=non_numeric_cols)
    # remake metadata pc data frame
    final_metadata_pc_data_frame_corrected = pd.concat([final_metadata_pc_data_frame[['SAMPLE']], numeric_vars, dummy_vars], axis=1)
    # calculate correlation matrix
    correlation_matrix = final_metadata_pc_data_frame_corrected.corr(numeric_only=True)
    # Print correlation matrix to standard output
    print(correlation_matrix)
    # make heatmap
    correlation_heatmap = sns.heatmap(correlation_matrix, cmap="YlGnBu", annot=False)
    # export heatmap to file
    plt.savefig(output_prefix + "_correlation_heatmap.png", dpi=300, bbox_inches='tight')
    # export correlation matrix to file
    correlation_matrix.to_csv(output_prefix + "_correlation_matrix.csv",index=True)
    # close figure
    plt.clf()

# main subroutine
def main():
    # Parse input and output arguments.
    args = parse_args()
    # import metadata data frame
    input_metadata_df=pd.read_csv(args.input_metadata)
    # import pc data frames
    input_pc_data_frame_list=[0] * len(args.input_pcs)
    for idx, i in args.input_pcs:
        input_pc_data_frame_list[idx]=pd.read_csv(i)
    # merge pc data frames with metadata data frame on SAMPLE column
    output_merged_metadata_pc_df=join_metadata_pc_tables(input_metadata_df,input_pc_data_frame_list,args.pc_counts)
    # output merged metadata/pc data frame to file
    output_merged_metadata_pc_df.to_csv(args.output_prefix + ".csv", index=False)
    # make covariate/PC correlation heatmap and save to file
    make_covariate_correlation_heatmap(output_merged_metadata_pc_df,args.output_prefix)
    
# run main subroutine
if __name__ == "__main__":
    main()
    