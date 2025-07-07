# script to make boxplots/violinplots of allele-specific QTL results (presence or absence of variant vs. phenotype)
# also generates output
# import modules first
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import time
import psutil
import argparse
import os

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Generate boxplot/violinplot phenotype distribution visualizations and merged tables for particular phenotype/variant combinations.")
    # import list of variant/phenotype combinations with tabs separating phenotype from CSV file
    parser.add_argument("--input_combinations", required=True, help="Headerless CSV file with list of comma-separated phenotype/variant combinations to be analyzed.")
    # import genetic data table
    parser.add_argument("--genetic_data", required=True, help="Path to genetic data input file; format described in CARDlongread_data_standardization repository.")
    # import methylation data table
    parser.add_argument("--methylation_data", required=True, help="Path to methylation data input file; format described in CARDlongread_data_standardization repository.")
    # set prefix for output files
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (plot and table per variant/phenotype combination).")
    # analyze unphased input data
    parser.add_argument("--unphased_input", action=argparse.BooleanOptionalAction, default=False, required=False, help="Visualize unphased QTL data outputs.")
    # compare phenotype for haplotype across all genotype haplotypes (trans comparison)
    parser.add_argument('--all_haps_comparison', action=argparse.BooleanOptionalAction, default=False, help="Compare all genetic against all methylation haplotypes (genetic H1/H2 against methylation H1/H2).")
    # add option for violinplot instead of boxplot
    parser.add_argument('--violin_plot', action=argparse.BooleanOptionalAction, default=False, dest="violin_plot", help="Show violin plots instead of box plots (optional; default false)")
    # add option for stripplot instead of swarmplot (in case of excessive data points)
    parser.add_argument('--strip_plot', action=argparse.BooleanOptionalAction, default=False, dest="strip_plot", help="Show strip plots instead of swarm plots inside box/violin plots (optional; default false)")
    return parser.parse_args()

# define phenotype plotting function
def plot_phenotypes(merged_data,variant_name,phenotype_name,violin_plot,strip_plot,output_prefix):
    # initialize figure and axes with matplotlib.plt
    fig, ax = plt.subplots()
    # different plots depending on whether violin_plot variable is set to true or not
    if violin_plot is True:
        ax = sb.violinplot(data=merged_data,x=variant_name,y=phenotype_name,color='white',inner="quartile")
    else:
        ax = sb.boxplot(data=merged_data,x=variant_name,y=phenotype_name,color='white')
    # different overlayed point plots depending on whether strip_plot is set 
    if strip_plot is True:
        ax = sb.stripplot(data=merged_data,x=variant_name,y=phenotype_name,ax=ax,dodge=True)
    else:
        ax = sb.swarmplot(data=merged_data,x=variant_name,y=phenotype_name,ax=ax,dodge=True)
    # print phenotype plot to png file
    fig.savefig(output_prefix + "_" + phenotype_name + "_" + variant_name + "_phenotype_plot.png", format='png', dpi=200, bbox_inches='tight')
    # close figure
    fig.clf()

# main script subroutine    
def main():
    # Parse input and output arguments.
    args = parse_args()
    
    # Start timing and memory tracking
    start_time = time.time()
    process = psutil.Process()
    
    # import input combinations
    input_combinations_df=pd.read_csv(args.input_combinations,header=None)
    # import genetic data
    genetic_data_df=pd.read_csv(args.genetic_data)
    # import methylation data
    methylation_data_df=pd.read_csv(args.methylation_data)
    
    # loop through input combinations by data frame row
    for index, row in input_combinations_df.iterrows():
        # define current phenotype based on current row (first column is phenotype)
        current_phenotype=row[0]
        # define current variant based on current row (second column is variant)
        current_variant=row[1]
        if (args.unphased_input is True):
            # subset genetic data on current variant of interest
            genetic_data_subset_df=genetic_data_df[['SAMPLE',current_variant]]
            # subset methylation data on current region of interest
            methylation_data_subset_df=methylation_data_df[['SAMPLE',current_phenotype]]
        else:
            # subset genetic data on current variant of interest
            genetic_data_subset_df=genetic_data_df[['SAMPLE','HAPLOTYPE',current_variant]]
            # subset methylation data on current region of interest
            methylation_data_subset_df=methylation_data_df[['SAMPLE','HAPLOTYPE',current_phenotype]]
        # cis or trans comparison depending on input setting
        if (args.all_haps_comparison is True) or (args.unphased_input is True):
            # merge genetic and methylation data on sample column alone
            merged_subset_df=pd.merge(genetic_data_subset_df,methylation_data_subset_df,on=['SAMPLE'])
        else:
            # merge genetic and methylation data on sample and haplotype columns
            merged_subset_df=pd.merge(genetic_data_subset_df,methylation_data_subset_df,on=['SAMPLE','HAPLOTYPE'])
        # make plot of current variant/phenotype pair
        # generate boxplots/violinplots with overlaid strip or swarmplots
        plot_phenotypes(merged_subset_df,current_variant,current_phenotype,args.violin_plot,args.strip_plot,args.output_prefix)
        # save merged data as csv
        merged_subset_df.to_csv(args.output_prefix + "_" + current_phenotype + "_" + current_variant + "_merged_table.csv")
    
    
if __name__ == "__main__":
    main()