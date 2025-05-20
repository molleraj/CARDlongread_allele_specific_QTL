#!/usr/bin/env python3
# script to generate PCs from different input data types for downstream analysis (e.g., genetic variants, regional methylation, gene expression, etc.)
# based on calculate_pcs.py by Spencer Grant

import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import argparse

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Perform exploratory data analysis (EDA) by running PCA for 20 PCs on input data (genetics, methylation, or expression) and generating a scree plot to assist selection of PCs.")
    # specify input data type (genetics, methylation, or expression)
    parser.add_argument("--input_type", required=True, help="Input file type (currently genetics, methylation, or expression). Genetics input is Plink eigenvalue/eigenvector files generated with the --pca 20 option, while methylation input is a methylation BED file and expression input is a normalized expression BED file.")
    # specify path to input file
    parser.add_argument("--input", required=True, help="Path to input genetics file prefix (for both Plink eigenvalue/eigenvector) or methylation/expression input file.")
    # specify output prefix
    parser.add_argument("--output_prefix", required=True, help="Prefix for PC and scree plot output files.")
    # specify PC prefix
    parser.add_argument("--pc_prefix", required=True, help="Prefix for PC names (e.g., METH_PC1).")
    # specify cutoff line for cumulative variance explained
    parser.add_argument("--cumulative_variance_explained_cutoff", required=False, default=0.70, help="Horizontal cutoff line for cumulative variance explained in scree plot (default 0.70). PC that exceeds threshold printed to standard output.")
    # specify missing information filter
    # specify plot title
    parser.add_argument("--plot_title", required=False, default="PCA scree plot", help="Title for output scree plot.")
    # add alternative sample name list
    parser.add_argument("--new_sample_names", required=False, default=None, help="Path to list of new sample names for output PC file. Must have same number of samples as input data file.")
    # return arguments
    return parser.parse_args()
    
# preprocess different data types
# preprocess genetics data input eigenvalue and eigenvector
def preprocess_genetic_inputs(plink_eigenvalues,plink_eigenvectors,pc_prefix):
    # convert plink eigenvalues into variance explained ratios
    plink_variance_explained_ratios = plink_eigenvalues/sum(plink_eigenvalues[0])
    # make copy of original eigenvectors DF
    plink_preprocessed_eigenvectors=plink_eigenvectors.copy()
    # rename plink eigenvector #IID or IID column to SAMPLE
    if '#IID' in plink_preprocessed_eigenvectors.columns:
        plink_preprocessed_eigenvectors.rename(columns={'#IID': 'SAMPLE'}, inplace=True)
    elif 'IID' in plink_preprocessed_eigenvectors.columns:
        plink_preprocessed_eigenvectors.rename(columns={'IID': 'SAMPLE'}, inplace=True)
    # get PC columns from eigenvectors column
    plink_eigenvector_PC_columns=plink_eigenvectors.columns[plink_eigenvectors.columns.str.contains('PC')]
    # relabel plink eigenvectors with PC prefix (e.g., GENETIC_)
    plink_eigenvectors_new_PC_names=[(i,pc_prefix+i) for i in plink_eigenvector_PC_columns.values]
    plink_preprocessed_eigenvectors.rename(columns=dict(plink_eigenvectors_new_PC_names), inplace=True)
    # return plink variance explained ratios and relabeled eigenvectors (PCs)
    return(plink_variance_explained_ratios,plink_preprocessed_eigenvectors)

# currently handles methylation and expression
# genetic PCs/eigenvalues come from Plink --pca
def preprocess_data_table(input_data_type,input_data_frame):
    # methylation data
    if (input_data_type == "methylation"):
        # handle methylation BED file
        # get columns with methylation levels - include avgMod or modFraction in name
        methylation_columns=input_data_frame.columns[input_data_frame.columns.str.contains('Mod|mod')]
        # make new data frame with methylation data only
        methylation_only_df=input_data_frame[methylation_columns]
        # correct names of methylation columns
        methylation_only_df.columns = methylation_only_df.columns.str.removeprefix("avgMod_")
        methylation_only_df.columns = methylation_only_df.columns.str.removesuffix("_GRCh38_1_modFraction") 
        # transpose to get samples as rows
        methylation_only_df = methylation_only_df.T
        # impute missing data as median - concern that mean maybe reflects outliers skewing distribution in one direction or other
        imputer = SimpleImputer(missing_values = np.nan, strategy ='median')
        # fit imputer to data and transform into final preprocessed data frame
        preprocessed_data_frame=pd.DataFrame(imputer.fit_transform(methylation_only_df))
        # relabel indices as samples
        preprocessed_data_frame.index=methylation_only_df.index
        # convert index to SAMPLE column
        # preprocessed_data_frame.reset_index(names="SAMPLE", inplace=True)
    # expression data
    elif (input_data_type == "expression"):
        # handle gene expression file
        # get columns with expression levels - ignore first four columns (chr, start, end, phenotype_id)
        expression_columns=input_data_frame.columns[4:]
        # different method to do this depending on whether TPM CSV input or extended BED file
        # make new data frame with expression data only
        expression_only_df=input_data_frame[expression_columns]
        # correct names of expression columns
        # expression_only_df.columns = expression_only_df.columns.str.removeprefix("avgMod_")
        # transpose to get samples as rows
        expression_only_df = expression_only_df.T
        # impute missing data as median - concern that mean maybe reflects outliers skewing distribution in one direction or other
        imputer = SimpleImputer(missing_values = np.nan, strategy ='median')
        # fit imputer to data and transform into final preprocessed data frame
        preprocessed_data_frame=pd.DataFrame(imputer.fit_transform(expression_only_df))
        # relabel indices as samples
        preprocessed_data_frame.index=expression_only_df.index
        # convert index to SAMPLE column
        # preprocessed_data_frame.reset_index(names="SAMPLE", inplace=True)
    # function returns preprocessed data frame for PCA
    return(preprocessed_data_frame)
        
# run pca
def run_pca(preprocessed_data_frame,pc_prefix):
    # calculate first 20 PCs
    # set up PCA
    pca = PCA(n_components=20, random_state=3)
    # make PCs with fit_transform method
    # exclude sample column
    pcs = pca.fit_transform(preprocessed_data_frame)
    # make PC data frame
    # rename with input PC prefix
    df_pcs = pd.DataFrame(pcs, index=preprocessed_data_frame.index, columns=[pc_prefix+f"PC{i}" for i in range(1, 21)])
    # convert index to SAMPLE column in final PC data frame for export
    df_pcs.reset_index(names="SAMPLE", inplace=True)
    # get variance explained ratios from pca model attribute
    variance_explained_ratios=pca.explained_variance_ratio_
    # return PC dataframe and variance explained ratios
    return(df_pcs,variance_explained_ratios)
    
# make scree plot of cumulative variance explained with user set horizontal cutoff (default 0.7 or 70% cumulative variance explained cutoff)
# also print first PC with which cumulative variance explained exceeds cutoff
def make_scree_plot(variance_explained_ratios,cumulative_variance_explained_cutoff,output_prefix,plot_title):
    # plot cumulative variance explained ratio vs. principal component number
    # calculate cumulative variance explained ratios
    cumulative_variance_explained_ratios=variance_explained_ratios.cumsum()
    # print PC upon which cumulative variance explained exceeds threshold IF any PCs exceed threshold
    if (len(cumulative_variance_explained_ratios[cumulative_variance_explained_ratios>cumulative_variance_explained_cutoff]) > 0):
        # print(cumulative_variance_explained_ratios)
        # pc_exceeding_threshold = np.nanargmin(cumulative_variance_explained_ratios[cumulative_variance_explained_ratios>cumulative_variance_explained_cutoff])+1
        # use np.argmax to find first true value instead
        pc_exceeding_threshold = np.argmax(cumulative_variance_explained_ratios>cumulative_variance_explained_cutoff)+1
        print("Cumulative variance exceeds " + str(cumulative_variance_explained_cutoff) + " threshold starting with PC" + pc_exceeding_threshold.astype(str))
    else:
        print("Cumulative variance below threshold for all 20 PCs.")
    # make plot_data dataframe
    # note that 20 PCs shown as 20 PCs determined in run_pca
    plot_data = pd.DataFrame({
        'PC': np.arange(1, 21),
        'CumulativeVarianceExplained': cumulative_variance_explained_ratios
    })
    # set figure size
    sns.set(rc={'figure.figsize':(12, 8)})
    # draw point plot with line segments joining points
    scree_plot = sns.pointplot(x='PC', y='CumulativeVarianceExplained', data=plot_data)
    # set plot style
    sns.set_style("whitegrid")
    # set plot context
    sns.set_context("talk")
    # turn off plot grid
    plt.grid(False)
    # define title based on plot_title input
    scree_plot.set_title(plot_title)
    # label y axis
    scree_plot.set_ylabel('Cumulative Variance Explained')
    # label x axis
    scree_plot.set_xlabel('Principal Components')
    # y axis limits of 0 and maximum cumulative variance explained plus 0.05 (5%)
    plt.ylim(0, cumulative_variance_explained_ratios.max() + 0.05)
    # plot horizontal line through variance explained cutoff
    plt.axhline(y=cumulative_variance_explained_cutoff,color="red")
    # save output figure
    plt.savefig(output_prefix + "_screeplot.png", dpi=300, bbox_inches='tight')
    # close figure
    plt.clf()
    
# main subroutine
def main():
    # Parse input and output arguments.
    args = parse_args()
    # import new sample names if specified
    if args.new_sample_names is not None:
        new_sample_names_df = pd.read_csv(args.new_sample_names, header=None)
    # import input data and handle accordingly
    if (args.input_type == 'genetics'):
        # load eigenvalues file
        genetics_eigenvalues = pd.read_csv(args.input + ".eigenval", header=None)
        # load eigenvectors (genetic PCs) file
        genetics_eigenvectors = pd.read_csv(args.input + ".eigenvec", sep="\t")
        # convert to data to relabeled PC table and variance explained ratios array
        (genetics_variance_explained_ratios,genetics_preprocessed_eigenvectors)=preprocess_genetic_inputs(genetics_eigenvalues,genetics_eigenvectors,args.pc_prefix)
        # generate scree plot
        make_scree_plot(genetics_variance_explained_ratios[0],args.cumulative_variance_explained_cutoff,args.output_prefix,args.plot_title)
        # relabel samples if specified
        if args.new_sample_names is not None:
            genetics_preprocessed_eigenvectors['SAMPLE']=new_sample_names_df[0]
        # export new PC table to csv; exclude indexes from CSV
        genetics_preprocessed_eigenvectors.to_csv(args.output_prefix + ".csv", index=False)
    elif (args.input_type == 'methylation'):
        # load methylation BED file
        methylation_input_df = pd.read_csv(args.input, sep="\t", na_values=['.'])
        # preprocess data
        methylation_preprocessed_data_frame = preprocess_data_table('methylation',methylation_input_df)
        # run PCA
        (methylation_df_pcs,methylation_variance_explained_ratios)=run_pca(methylation_preprocessed_data_frame,pc_prefix)
        # generate scree plot
        make_scree_plot(methylation_variance_explained_ratios,args.cumulative_variance_explained_cutoff,args.output_prefix,args.plot_title)
        # relabel samples if specified
        if args.new_sample_names is not None:
            methylation_df_pcs['SAMPLE']=new_sample_names_df[0]
        # export new PC table to csv; exclude indexes from CSV
        methylation_df_pcs.to_csv(args.output_prefix + ".csv",index=False)
    elif (args.input_type == 'expression'):
        # load normalized expression BED file
        expression_input_df = pd.read_csv(args.input, sep="\t", na_values=['.'])
        # preprocess data
        expression_preprocessed_data_frame = preprocess_data_table('expression',expression_input_df)
        # run PCA
        (expression_df_pcs,expression_variance_explained_ratios)=run_pca(expression_preprocessed_data_frame,pc_prefix)
        # generate scree plot
        make_scree_plot(expression_variance_explained_ratios,args.cumulative_variance_explained_cutoff,args.output_prefix,args.plot_title)
        # relabel samples if specified
        if args.new_sample_names is not None:
            expression_df_pcs['SAMPLE']=new_sample_names_df[0]
        # export new PC table to csv; exclude indexes from CSV
        expression_df_pcs.to_csv(args.output_prefix + ".csv", index=False)
# run main subroutine
if __name__ == "__main__":
    main()