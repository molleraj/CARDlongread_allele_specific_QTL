# script to join my results with Ken's results
# pull top (most significant p value) hits per region and compare
# import modules first
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import time
import psutil
import argparse
import os

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Compare allele-specific QTL output to tensorQTL output on common methylation region. Print common regions and those unique to each.")
    parser.add_argument("--tensor_QTL", required=True, help="Path to input tensorQTL tsv file for comparisons; 'phenotype_id	variant_id	af	pval_nominal	slope	slope_se	pval_perm	bh_fdr	qval	Chromosome	TOP SV ID	TOP SV Causal Post Probablity	TOP SNV ID	TOP SNV Causal Post Probablity' as header.")
    parser.add_argument("--allele_specific_QTL", required=True, help="Path to input allele-specific QTL file for comparisons; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.")
    parser.add_argument("--region_type", required=True, help="Region type (CGI for CpG islands, GB for gene bodies, PROM for promoters).")
    parser.add_argument("--plot_title", required=False, help="Title for beta comparison histogram and scatterplot.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files with tensorQTL only, allele-specific QTL only, and common hits.")
    return parser.parse_args()

# subroutine to merge tables
def merge_tables(tensorQTL_df,allele_spec_QTL_df,region_type):
    if (region_type == "PROM"):
        # if region type is PROM, then just use name of promoter
        # only split into three additional elements - leave underscore in promoter name alone
        tensorQTL_df['pheno_name']=tensorQTL_df['phenotype_id'].str.split("_",n=3,expand=True)[3]
        # promoter name is just outcome variable
        allele_spec_QTL_df['pheno_name']=allele_spec_QTL_df['outcome']
        # merge on pheno_name
        merged_QTL_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name')
        # get tensorQTL or allele-specific QTL only hits
        # use outer join
        tensor_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name',how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name',how="outer",indicator=True).query('_merge=="right_only"')
    else:
        # split tensorQTL phenotype_id by underscore and convert respective fields into chromosome, start, and end for comparison
        tensorQTL_df['pheno_chrom']=tensorQTL_df['phenotype_id'].str.split("_",expand=True)[0]
        tensorQTL_df['pheno_start']=tensorQTL_df['phenotype_id'].str.split("_",expand=True)[1]
        tensorQTL_df['pheno_end']=tensorQTL_df['phenotype_id'].str.split("_",expand=True)[2]
        # do same split for allele-specific QTL hits based on region name (outcome column)
        allele_spec_QTL_df['pheno_chrom']=allele_spec_QTL_df['outcome'].str.split("_",expand=True)[1]
        allele_spec_QTL_df['pheno_start']=allele_spec_QTL_df['outcome'].str.split("_",expand=True)[2]
        allele_spec_QTL_df['pheno_end']=allele_spec_QTL_df['outcome'].str.split("_",expand=True)[3]
        # now merge tables on above three fields
        merged_QTL_df=tensorQTL_df.merge(allele_spec_QTL_df,on=['pheno_chrom','pheno_start','pheno_end'])
        # get tensorQTL or allele-specific QTL only hits
        # use outer join
        tensor_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on=['pheno_chrom','pheno_start','pheno_end'],how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on=['pheno_chrom','pheno_start','pheno_end'],how="outer",indicator=True).query('_merge=="right_only"')
    # compare beta sign for matches
    merged_QTL_df['beta_signs_match']=(np.sign(merged_QTL_df['slope'])==np.sign(merged_QTL_df['beta']))
    # print number of concordant signs, discordant signs, and total samples
    print("Concordant beta signs: ",merged_QTL_df['beta_signs_match'].value_counts()[0])
    print("Discordant beta signs: ",merged_QTL_df['beta_signs_match'].value_counts()[1])
    print("Total common regions: ",len(merged_QTL_df))
    return(merged_QTL_df,tensor_QTL_only_df,allele_specific_QTL_only_df)
    
# get hits only in one or other data frame    
# def get_excluded_hits():
def common_hits_visualizations(common_QTL_hits_df,output_prefix,plot_title):
    # output beta histogram
    fig, ax = plt.subplots()
    # overlay unphased betas over phased betas
    # phased betas
    ax = sb.histplot(data=common_QTL_hits_df,x=beta,label="Phased")
    # unphased betas
    sb.histplot(data=common_QTL_hits_df,x=slope,label="Unphased",ax=ax)
    # set axis labels
    ax.set(xlabel="Beta",ylabel="Frequency")
    # include legend
    plt.legend()
    # save histogram figure
    fig.savefig(output_prefix + "_beta_histogram.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    # output unphased vs. phased betas as scatterplot with linear regression line
    # print scatterplot with regression line
    fig, ax = plt.subplots()
    ax = sb.regplot(data=common_QTL_hits_df,x=slope,y=beta)
    # set axis labels
    ax.set(xlabel="Unphased beta",ylabel="Phased beta")
    if plot_title is not None:
        ax.set(title=plot_title)
    # save scatterplot figure
    fig.savefig(output_prefix + "_beta_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    
# main script subroutine    
def main():
    # Parse input and output arguments.
    args = parse_args()
    
    # Start timing and memory tracking
    start_time = time.time()
    process = psutil.Process()
    
    # Load input files
    # note that tensorQTL output is tab separated
    tensorQTL_df = pd.read_csv(args.tensor_QTL,sep="\t")
    allele_spec_QTL_df = pd.read_csv(args.allele_specific_QTL)
    
    # Merge input files
    (merged_QTL_df,tensor_QTL_only_df,allele_specific_QTL_only_df)=merge_tables(tensorQTL_df,allele_spec_QTL_df,args.region_type)
    
    # output separate data frames
    merged_QTL_df.to_csv(args.output_prefix + "_common.csv",index=False)
    tensor_QTL_only_df.to_csv(args.output_prefix + "_tensor_QTL_only.csv",index=False)
    allele_specific_QTL_only_df.to_csv(args.output_prefix + "_allele_specific_QTL_only.csv",index=False)
    
    # Print runtime and max RAM usage
    end_time = time.time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    print(f"Max RAM Usage: {process.memory_info().rss / (1024 ** 2):.2f} MB")

if __name__ == "__main__":
    main()