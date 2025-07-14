# script to filter QTL output by p-value using BH FDR correction

import pandas as pd
import numpy as np
import time
import psutil
from statsmodels.stats import multitest
import argparse
import os

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Perform multiple hypothesis testing false-discovery rate (FDR) correction for allele-specific methylation QTL output from all chromosomes.")
    parser.add_argument("--input", required=True, help="Path to input QTL file for all chromosomes; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq,intercept,overall_MAF,H1_MAF,H2_MAF,H1_missing_meth,H2_missing_meth' as header.")
    parser.add_argument("--output", required=True, help="Path to the output QTL file with additional two fields - rejected and pvalue-corrected.")
    parser.add_argument("--rejected", action=argparse.BooleanOptionalAction, default=False, required=False, help="Set this option to only print results where the null hypothesis was rejected (i.e., statistically significant associations).")
    parser.add_argument("--drop_na", action=argparse.BooleanOptionalAction, default=False, required=False, help="Drop any cases where p value is listed as NA.")
    parser.add_argument("--top_hit_per_gene", action=argparse.BooleanOptionalAction, default=False, required=False, help="Print only most significant hit per gene.")
    parser.add_argument("--top_hit_per_region", action=argparse.BooleanOptionalAction, default=False, required=False, help="Print only most significant hit per methylation region.")
    # add argument to choose top by p value or R2
    parser.add_argument("--top_by_property",choices=["p_value", "r2"],default="p_value",required=False,help="Choose top hit by this property (p value or R2 value).")
    parser.add_argument("--unique_hits", action=argparse.BooleanOptionalAction, default=False, required=False, help="Print only unique set of included tests.")
    return parser.parse_args()

# get most significant QTL hit per gene    
def get_top_hit_per_gene(QTL_df,top_by_property):
    # get unique gene list
    unique_gene_list=np.unique(QTL_df['gene'].astype(str))
    # initialize output data frame
    top_hit_per_gene_df = pd.DataFrame({'gene' : unique_gene_list},columns=QTL_df.columns)
    # for each gene, select most significant hit
    for idx, i in enumerate(unique_gene_list):
        # get QTL data frame line with smallest p value for given gene
        QTL_df_current_gene=QTL_df[QTL_df['gene']==i]
        # get first line if multiple hits
        top_hit_per_current_gene=QTL_df_current_gene[QTL_df_current_gene[top_by_property] == min(QTL_df_current_gene[top_by_property])].iloc[0,:]
        # set current row of top_hit_per_gene_df to first hit line per gene (amongst most significant hits)
        top_hit_per_gene_df.iloc[idx,:]=top_hit_per_current_gene
    # return new data frame only with most significant hit for each gene
    return(top_hit_per_gene_df)
    
# get most significant QTL hit per region    
def get_top_hit_per_region(QTL_df,top_by_property):
    # get unique region list
    unique_region_list=np.unique(QTL_df['outcome'].astype(str))
    # initialize output data frame
    top_hit_per_region_df = pd.DataFrame({'outcome' : unique_region_list},columns=QTL_df.columns)
    # for each region, select most significant hit
    for idx, i in enumerate(unique_region_list):
        # get QTL data frame line with smallest p value for given region
        QTL_df_current_region=QTL_df[QTL_df['outcome']==i]
        # get first line if multiple hits
        top_hit_per_current_region=QTL_df_current_region[QTL_df_current_region[top_by_property] == min(QTL_df_current_region[top_by_property])].iloc[0,:]
        # set current row of top_hit_per_region_df to first hit line per region (amongst most significant hits)
        top_hit_per_region_df.iloc[idx,:]=top_hit_per_current_region
    # return new data frame only with most significant hit for each region
    return(top_hit_per_region_df)

# get unique set of QTL hits
def get_unique_QTL_hits(QTL_df):
    # get unique rows, excluding gene
    unique_QTL_hits_df=QTL_df.loc[QTL_df.drop(columns='gene').drop_duplicates().index] 
    # return data frame only with unique pairs of outcomes (methylation) and predictors (genotypes)
    return(unique_QTL_hits_df)

# main script subroutine    
def main():
    # Parse input and output arguments.
    args = parse_args()
    
    # Start timing and memory tracking
    start_time = time.time()
    process = psutil.Process()
    
    # Load input file
    input_QTL_df = pd.read_csv(args.input)
    
    # print total tests considered and file name
    print("Input QTL file is ",args.input)
    print("Total tests processed: ",len(input_QTL_df))
    
    # drop na p values if drop_na selected
    # if na p values are kept, all corrected p-values are na
    if args.drop_na is True:
        input_QTL_df = input_QTL_df[input_QTL_df.p_value.isna() != True]
    
    # perform Benjamini-Hochberg multiple hypothesis testing correction
    # BH correction is default. Also unsorted p-values is default
    # fdr correction output is two arrays - first rejection (True or False for null hypothesis) and second corrected p-value
    input_QTL_df['rejected']=multitest.fdrcorrection(input_QTL_df.p_value)[0]
    input_QTL_df['pvalue-corrected']=multitest.fdrcorrection(input_QTL_df.p_value)[1]
    
    # sequentially filter input_QTL_df first for unique hits, then h0 rejected (significant) hits, then top hit per region, and then top hit per gene
    # note cascading filtering effects if all selected
    
    if args.unique_hits is True:
        # filter only for unique hits
        # may remove gene names (only keeps first of hit found multiple times)
        input_QTL_df=get_unique_QTL_hits(input_QTL_df)
        print("Total unique tests: ", len(input_QTL_df))
    
    if args.rejected is True:
        # only print results where null hypothesis is rejected
        input_QTL_df=input_QTL_df[input_QTL_df['rejected']==True]
        print("Total significant tests: ", len(input_QTL_df))
        
    if args.top_hit_per_region is True:
        # filter for only top hit for each gene
        input_QTL_df=get_top_hit_per_region(input_QTL_df,args.top_by_property)
        print("Total significant regions: ",len(input_QTL_df))
    
    if args.top_hit_per_gene is True:
        # filter for only top hit for each gene
        input_QTL_df=get_top_hit_per_gene(input_QTL_df,args.top_by_property)
        print("Total significant genes of interest: ",len(input_QTL_df))
    
    # save output file; no indexes
    print("Output QTL file is ",args.output)
    input_QTL_df.to_csv(args.output,index=False)
    
    # Print runtime and max RAM usage
    end_time = time.time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    print(f"Max RAM Usage: {process.memory_info().rss / (1024 ** 2):.2f} MB")
    
if __name__ == "__main__":
    main()