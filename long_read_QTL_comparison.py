# script to join my results with Ken's tensorQTL results or those of an additional allele-specific run
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
    parser = argparse.ArgumentParser(description="Compare allele-specific QTL output to tensorQTL or additional allele-specific QTL output on common methylation region. Print common regions and those unique to each.")
    # argument for whether to compare against another allele specific QTL run or a tensor QTL run
    parser.add_argument("--comparison_type",choices=["allele_specific_QTL_unphased", "allele_specific_QTL_phased","allele_specific_QTL_dephased","tensor_QTL"],default="tensor_QTL",required=False,help="Comparison against tensorQTL hits or additional allele-specific QTL run.")
    # argument for how to join allele specific QTL tables (e.g., for significant unique vs. all unique hits between runs)
    parser.add_argument("--allele_specific_join_on",choices=["region","variant","both"],default="region",required=False,help="Fields on which to join allele-specific QTL tables - methylation region (outcome), genetic variant (predictor), or both.")
    parser.add_argument("--tensor_QTL", required=False, help="Path to input tensorQTL tsv file for comparisons; 'phenotype_id	variant_id	af	pval_nominal	slope	slope_se	pval_perm	bh_fdr	qval	Chromosome	TOP SV ID	TOP SV Causal Post Probablity	TOP SNV ID	TOP SNV Causal Post Probablity' as header.")
    # require at least one argument for allele specific QTL input;
    parser.add_argument("--allele_specific_QTL", required=True, nargs="+", help="Path to input allele-specific QTL file(s) for comparisons; either single QTL input or forced unphased/phased in that order; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as minimal header.")
    parser.add_argument("--region_type", required=True, help="Region type (CGI for CpG islands, GB for gene bodies, PROM for promoters).")
    parser.add_argument("--plot_title", required=False, help="Title for beta comparison histogram and scatterplot.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files with tensorQTL only, allele-specific QTL only, and common hits.")
    # add argument for methylation regions to include
    parser.add_argument("--include_regions", required=False, help="BED file containing list of regions to include in comparison or list of promoter names in the case of promoters (filter out all others).")
    return parser.parse_args()

# subroutine to merge tensorQTL with allele-specific QTL tables
def merge_tensor_allele_spec_tables(tensorQTL_df,allele_spec_QTL_df,include_regions,region_type):
    if (region_type == "PROM"):
        # if region type is PROM, then just use name of promoter
        # only split into three additional elements - leave underscore in promoter name alone
        tensorQTL_df['pheno_name']=tensorQTL_df['phenotype_id'].str.split("_",n=3,expand=True)[3]
        # promoter name is just outcome variable
        allele_spec_QTL_df['pheno_name']=allele_spec_QTL_df['outcome']
        # filter each on include regions
        if (include_regions is not None):
            include_regions.columns=['pheno_name']
            tensorQTL_df=tensorQTL_df.merge(include_regions,on='pheno_name',how="inner")
            allele_spec_QTL_df=allele_spec_QTL_df.merge(include_regions,on='pheno_name',how="inner")
        # merge on pheno_name
        merged_QTL_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name')
        # get tensorQTL or allele-specific QTL only hits
        # use outer join
        tensor_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name',how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL_only_df=tensorQTL_df.merge(allele_spec_QTL_df,on='pheno_name',how="outer",indicator=True).query('_merge=="right_only"')
    else:
        # split tensorQTL phenotype_id by underscore and convert respective fields into chromosome, start, and end for comparison
        tensorQTL_df['pheno_chrom']=tensorQTL_df['phenotype_id'].str.split("_",expand=True)[0]
        tensorQTL_df['pheno_start']=pd.to_numeric(tensorQTL_df['phenotype_id'].str.split("_",expand=True)[1])
        tensorQTL_df['pheno_end']=pd.to_numeric(tensorQTL_df['phenotype_id'].str.split("_",expand=True)[2])
        # do same split for allele-specific QTL hits based on region name (outcome column)
        allele_spec_QTL_df['pheno_chrom']=allele_spec_QTL_df['outcome'].str.split("_",expand=True)[1]
        allele_spec_QTL_df['pheno_start']=pd.to_numeric(allele_spec_QTL_df['outcome'].str.split("_",expand=True)[2])
        allele_spec_QTL_df['pheno_end']=pd.to_numeric(allele_spec_QTL_df['outcome'].str.split("_",expand=True)[3])
        # filter each on include regions
        if (include_regions is not None):
            include_regions.columns=['pheno_chrom','pheno_start','pheno_end']
            tensorQTL_df=tensorQTL_df.merge(include_regions,on=['pheno_chrom','pheno_start','pheno_end'],how="inner")
            allele_spec_QTL_df=allele_spec_QTL_df.merge(include_regions,on=['pheno_chrom','pheno_start','pheno_end'],how="inner")
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
    if (len(merged_QTL_df['beta_signs_match'].value_counts()) > 1):
        print("Discordant beta signs: ",merged_QTL_df['beta_signs_match'].value_counts()[1])
    print("Total common regions: ",len(merged_QTL_df))
    return(merged_QTL_df,tensor_QTL_only_df,allele_specific_QTL_only_df)

# subroutine to merge two allele-specific QTL tables
def merge_two_allele_spec_tables(allele_spec_QTL_df1,allele_spec_QTL_df2,include_regions,region_type,join_type):
    if (region_type == "PROM"):
        if (include_regions is not None):
            include_regions.columns=['pheno_name']
            # promoter name is just outcome variable
            allele_spec_QTL_df1['pheno_name']=allele_spec_QTL_df1['outcome']
            allele_spec_QTL_df2['pheno_name']=allele_spec_QTL_df2['outcome']
            allele_spec_QTL_df1=allele_spec_QTL_df1.merge(include_regions,on='pheno_name',how="inner")
            allele_spec_QTL_df2=allele_spec_QTL_df2.merge(include_regions,on='pheno_name',how="inner")
    else:
        if (include_regions is not None):
            include_regions.columns=['pheno_chrom','pheno_start','pheno_end']
            # do same split for allele-specific QTL 1 hits based on region name (outcome column)
            allele_spec_QTL_df1['pheno_chrom']=allele_spec_QTL_df1['outcome'].str.split("_",expand=True)[1]
            allele_spec_QTL_df1['pheno_start']=pd.to_numeric(allele_spec_QTL_df1['outcome'].str.split("_",expand=True)[2])
            allele_spec_QTL_df1['pheno_end']=pd.to_numeric(allele_spec_QTL_df1['outcome'].str.split("_",expand=True)[3])
            # do same split for allele-specific QTL 2 hits based on region name (outcome column)
            allele_spec_QTL_df2['pheno_chrom']=allele_spec_QTL_df2['outcome'].str.split("_",expand=True)[1]
            allele_spec_QTL_df2['pheno_start']=pd.to_numeric(allele_spec_QTL_df2['outcome'].str.split("_",expand=True)[2])
            allele_spec_QTL_df2['pheno_end']=pd.to_numeric(allele_spec_QTL_df2['outcome'].str.split("_",expand=True)[3])
            allele_spec_QTL_df1=allele_spec_QTL_df1.merge(include_regions,on=['pheno_chrom','pheno_start','pheno_end'],how="inner")
            allele_spec_QTL_df2=allele_spec_QTL_df2.merge(include_regions,on=['pheno_chrom','pheno_start','pheno_end'],how="inner")
    # different joins depending on join_type
    if (join_type == "region"):
        # merge allele-specific QTL tables on outcome
        merged_QTL_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome'])
        # get allele-specific QTL run 1 or 2 only hits
        allele_specific_QTL1_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome'],how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL2_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome'],how="outer",indicator=True).query('_merge=="right_only"')
    elif (join_type == "variant"):
        # merge allele-specific QTL tables on outcome
        merged_QTL_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['predictor'])
        # get allele-specific QTL run 1 or 2 only hits
        allele_specific_QTL1_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['predictor'],how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL2_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['predictor'],how="outer",indicator=True).query('_merge=="right_only"')
    elif (join_type == "both"):
        # merge allele-specific QTL tables on outcome
        merged_QTL_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome','predictor'])
        # get allele-specific QTL run 1 or 2 only hits
        allele_specific_QTL1_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome','predictor'],how="outer",indicator=True).query('_merge=="left_only"')
        allele_specific_QTL2_only_df=allele_spec_QTL_df1.merge(allele_spec_QTL_df2,on=['outcome','predictor'],how="outer",indicator=True).query('_merge=="right_only"')
    # compare beta sign for matches
    merged_QTL_df['beta_signs_match']=(np.sign(merged_QTL_df['beta_x'])==np.sign(merged_QTL_df['beta_y']))
    # print number of concordant signs, discordant signs, and total samples
    print("Concordant beta signs: ",merged_QTL_df['beta_signs_match'].value_counts()[0])
    if (len(merged_QTL_df['beta_signs_match'].value_counts()) > 1):
        print("Discordant beta signs: ",merged_QTL_df['beta_signs_match'].value_counts()[1])
    print("Total common regions: ",len(merged_QTL_df))
    return(merged_QTL_df,allele_specific_QTL1_only_df,allele_specific_QTL2_only_df)
    
# get hits only in one or other data frame    
# def get_excluded_hits():
# make common hits beta histograms/scatterplot for comparison
def common_hits_visualizations(common_QTL_hits_df,comparison_type,output_prefix,plot_title):
    # set purple colors for plotting
    pcolors = sb.color_palette("Purples", 2)[::-1]
    # output beta histogram
    fig, ax = plt.subplots()
    # show unphased betas grouped side by side with phased betas
    # first create unphased and phased variables
    if (comparison_type == "tensor_QTL"):
        common_QTL_hits_df['Unphased']=common_QTL_hits_df['slope']
        common_QTL_hits_df['Phased']=common_QTL_hits_df['beta']
        # now melt on these variables
        common_QTL_hits_df_melted=common_QTL_hits_df.melt(value_vars=['Unphased', 'Phased'],
                    var_name='QTL type', 
                    value_name='Beta')
    elif (comparison_type == "allele_specific_QTL_unphased"):
        common_QTL_hits_df['Unphased']=common_QTL_hits_df['beta_x']
        common_QTL_hits_df['Phased']=common_QTL_hits_df['beta_y']
        # now melt on these variables
        common_QTL_hits_df_melted=common_QTL_hits_df.melt(value_vars=['Unphased', 'Phased'],
                    var_name='QTL type', 
                    value_name='Beta')
    elif (comparison_type == "allele_specific_QTL_phased"):
        common_QTL_hits_df['Phased run 1']=common_QTL_hits_df['beta_x']
        common_QTL_hits_df['Phased run 2']=common_QTL_hits_df['beta_y']
        # now melt on these variables
        common_QTL_hits_df_melted=common_QTL_hits_df.melt(value_vars=['Phased run 1', 'Phased run 2'],
                    var_name='QTL type', 
                    value_name='Beta')
    elif (comparison_type == "allele_specific_QTL_dephased"):
        common_QTL_hits_df['Dephased']=common_QTL_hits_df['beta_x']
        common_QTL_hits_df['Phased']=common_QTL_hits_df['beta_y']
        # now melt on these variables
        common_QTL_hits_df_melted=common_QTL_hits_df.melt(value_vars=['Dephased', 'Phased'],
                    var_name='QTL type', 
                    value_name='Beta')
    # then make grouped 
    sb.histplot(data=common_QTL_hits_df_melted,x='Beta',hue='QTL type',multiple="dodge",palette=pcolors)
    # set axis labels
    ax.set(xlabel="Beta",ylabel="Frequency")
    # set plot title if defined
    if plot_title is not None:
        ax.set(title=plot_title)
    # save histogram figure
    fig.savefig(output_prefix + "_beta_histogram.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    # output unphased vs. phased betas as scatterplot with linear regression line
    # print scatterplot with regression line
    fig, ax = plt.subplots()
    # different x and y depending on whether tensor QTL or another allele-specific QTL considered
    if (comparison_type == "tensor_QTL"):
        ax = sb.regplot(data=common_QTL_hits_df,x='slope',y='beta')
        # set axis labels
        ax.set(xlabel="Unphased beta",ylabel="Phased beta")
    elif (comparison_type == "allele_specific_QTL_unphased"):
        ax = sb.regplot(data=common_QTL_hits_df,x='beta_x',y='beta_y')
        # set axis labels
        ax.set(xlabel="Unphased beta",ylabel="Phased beta")
    elif (comparison_type == "allele_specific_QTL_phased"):
        ax = sb.regplot(data=common_QTL_hits_df,x='beta_x',y='beta_y')
        # set axis labels
        ax.set(xlabel="Phased run 1 beta",ylabel="Phased run 2 beta")
    elif (comparison_type == "allele_specific_QTL_dephased"):
        ax = sb.regplot(data=common_QTL_hits_df,x='beta_x',y='beta_y')
        # set axis labels
        ax.set(xlabel="Dephased beta",ylabel="Phased beta")
    # set plot title if defined
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
    # Load region filtering list
    if (args.include_regions is not None):
        regions_to_include_df=pd.read_csv(args.include_regions,sep="\t",header=None)
    else:
        regions_to_include_df=None
    # note that tensorQTL output is tab separated
    # only import if comparison type is "tensor_QTL"
    if (args.comparison_type == "tensor_QTL"):
        tensorQTL_df = pd.read_csv(args.tensor_QTL,sep="\t")
    if len(args.allele_specific_QTL) > 1:
        # forced unphased, unphased, or dephased
        allele_spec_QTL_df1 = pd.read_csv(args.allele_specific_QTL[0])
        # phased
        allele_spec_QTL_df2 = pd.read_csv(args.allele_specific_QTL[1])
    elif len(args.allele_specific_QTL) == 1:
        allele_spec_QTL_df = pd.read_csv(args.allele_specific_QTL[0])
    
    # Merge input files
    # depends on comparison type
    if (args.comparison_type == "tensor_QTL"):
        (merged_QTL_df,tensor_QTL_only_df,allele_specific_QTL_only_df)=merge_tensor_allele_spec_tables(tensorQTL_df,allele_spec_QTL_df,regions_to_include_df,args.region_type)
        # output separate data frames
        merged_QTL_df.to_csv(args.output_prefix + "_common.csv",index=False)
        tensor_QTL_only_df.to_csv(args.output_prefix + "_tensor_QTL_only.csv",index=False)
        allele_specific_QTL_only_df.to_csv(args.output_prefix + "_allele_specific_QTL_only.csv",index=False)
        # print total tensor QTL hits
        # print total allele specific phased QTL hits
        # print common QTL hits and percentage phased QTL hits
    else:
    # elif (args.comparison_type == "allele_specific_QTL"):
        (merged_QTL_df,allele_specific_QTL1_only_df,allele_specific_QTL2_only_df)=merge_two_allele_spec_tables(allele_spec_QTL_df1,allele_spec_QTL_df2,regions_to_include_df,args.region_type,args.allele_specific_join_on)
        # output separate data frames
        merged_QTL_df.to_csv(args.output_prefix + "_common.csv",index=False)
        allele_specific_QTL1_only_df.to_csv(args.output_prefix + "_allele_specific_QTL1_only.csv",index=False)
        allele_specific_QTL2_only_df.to_csv(args.output_prefix + "_allele_specific_QTL2_only.csv",index=False)
        # print total allele specific phased QTL run 1 hits
        # print total allele specific phased QTL run 2 hits
        # print common QTL hits and percentage phased QTL run 2 hits
    
    # make histogram and scatterplot of common hits' betas
    common_hits_visualizations(merged_QTL_df,args.comparison_type,args.output_prefix,args.plot_title)
    
    # Print runtime and max RAM usage
    end_time = time.time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    print(f"Max RAM Usage: {process.memory_info().rss / (1024 ** 2):.2f} MB")

if __name__ == "__main__":
    main()