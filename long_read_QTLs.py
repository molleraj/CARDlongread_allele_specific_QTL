import pandas as pd
import numpy as np
import time
import psutil
import statsmodels.api as sm
import statsmodels.formula.api as smf
import argparse
import os
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Perform QTL analysis with user-specified parameters. Input data formats described in CARDlongread_data_standardization repository.")
    parser.add_argument("--chromosome", required=True, help="Specify the chromosome to subset the ROI map (e.g., 'chr1').")
    parser.add_argument("--roi_map", required=True, help="Path to the ROI map file.")
    parser.add_argument("--genetic_data", required=True, help="Path to the genetic data file (e.g., '/path/to/genetic_data').")
    parser.add_argument("--methylation_data", required=True, help="Path to the methylation data file.")
    parser.add_argument("--genetic_map", required=True, help="Path to the genetic map file.")
    parser.add_argument("--methylation_map", required=True, help="Path to the methylation map file.")
    parser.add_argument("--metadata", required=True, help="Path to the metadata file.")
    parser.add_argument("--output", required=True, help="Output file destination for results.")
    parser.add_argument("--simulate_unphased", action=argparse.BooleanOptionalAction, default=False, required=False, help="Convert phased genetics/methylation data to unphased where possible (i.e., exclude NAs) and run QTL on unphased data.")
    parser.add_argument("--window_size", type=int, default=500000, help="Window size for defining gene regions (default: 500000).")
    # parser.add_argument("--minimum_genotypes",type=int, default=3, help="Minimum number of variant genotypes to run regression (default: 3).")
    parser.add_argument("--per_haplotype_missing_methylation_rate", default=0.95, help="Proportion of methylation calls missing per haplotype to consider for corresponding MAF exclusion (default: 0.95).")
    parser.add_argument("--per_haplotype_MAF", default=0.05, help="Minimum minor allele frequency per haplotype to run regression in the case of one haplotype with data (default: 0.05).")
    parser.add_argument("--overall_MAF", default=0.05, help="Minimum minor allele frequency for both haplotypes to run regression (default: 0.05).")
    return parser.parse_args()
    
def collapse_haps(genetic_data,methylation_data):
    # make collapsed haps genetic and methylation data frames
    collapsed_haps_genetic_data=pd.DataFrame(data={'SAMPLE' : np.unique(genetic_data['SAMPLE'])},columns=genetic_data.columns.drop('HAPLOTYPE'))
    collapsed_haps_methylation_data=pd.DataFrame(data={'SAMPLE' : np.unique(methylation_data['SAMPLE'])},columns=methylation_data.columns.drop('HAPLOTYPE'))
    # sum 0/1 across haps per variant, per sample
    # except NA cases (set genotype sum to be NA)
    for idx, i in enumerate(np.unique(genetic_data['SAMPLE'])):
        # set min_count to 2 to exclude NAs
        collapsed_haps_genetic_data.iloc[idx,1:]=genetic_data[genetic_data['SAMPLE']==i].drop(columns=['SAMPLE','HAPLOTYPE']).sum(min_count=2).values
    # average methylation data across haps per region, per sample
    # except NAs (set average to NA)
    for idx, i in enumerate(np.unique(methylation_data['SAMPLE'])):
        # set min_count to 2 to exclude NAs
        collapsed_haps_methylation_data.iloc[idx,1:]=methylation_data[methylation_data['SAMPLE']==i].drop(columns=['SAMPLE','HAPLOTYPE']).sum(min_count=2).values/2
    # return collapsed haps genetic and methylation data
    return(collapsed_haps_genetic_data,collapsed_haps_methylation_data)

def main():
    args = parse_args()
    
    # Start timing and memory tracking
    start_time = time.time()
    process = psutil.Process()

    # Load input files
    # roi map is a tab separated bed file
    roi_map = pd.read_csv(args.roi_map,sep="\t")
    genetic_map = pd.read_csv(args.genetic_map)
    methylation_map = pd.read_csv(args.methylation_map)
    metadata = pd.read_csv(args.metadata)
    # set types to reduce size of genetic/methylation data frames
    # set presence/absence as 8 bit integers allowing NAs (Int8)
    # genetic_data_types = defaultdict(lambda: "Int8", SAMPLE="str", HAPLOTYPE="str")
    # import genetic and methylation data given NA as na value
    genetic_data = pd.read_csv(args.genetic_data, na_values="NA") #(dtype=genetic_data_types)
    # methylation data has floats so don't use data_types
    methylation_data = pd.read_csv(args.methylation_data, na_values="NA")
    # collapse haps for genetic and methylation data if simulate_unphased set
    (unphased_genetic_data,unphased_methylation_data)=collapse_haps(genetic_data,methylation_data)
    
    # Subset ROI map by chromosome
    roi_map = roi_map[roi_map['CHROM'] == args.chromosome]
    results = []
    
    for _, gene_info in roi_map.iterrows():
        gene_of_interest = gene_info['NAME']
        chrom = gene_info['CHROM']
        start = gene_info['START'] - args.window_size
        end = gene_info['END'] + args.window_size

        # Get variants and probes in the defined region
        variants = genetic_map[(genetic_map['CHROM'] == chrom) & (genetic_map['START'].between(start, end))]
        probes = methylation_map[(methylation_map['CHROM'] == chrom) & ((methylation_map['START'] <= end) & (methylation_map['STOP'] >= start))]
        # set variant column to variant name
        # variants['VARIANT']=variants['NAME']
        variants=variants.rename(columns={'NAME':'VARIANT'})
        # set target column to probe name
        # probes['TARGET']=probes['NAME']
        probes=probes.rename(columns={'NAME':'TARGET'})
        # Load genetic and methylation data
        # genetic_file = os.path.join(args.genetic_data_dir, f"{chrom}.tsv")
        # if not os.path.exists(genetic_file):
            # print(f"Warning: Genetic file {genetic_file} not found. Skipping...")
            # continue
        
        # Process metadata (convert only non-numeric columns to dummy variables)
        covariates = metadata.drop(columns=['SAMPLE'])
        # for debugging
        # print(covariates)
        non_numeric_cols = covariates.select_dtypes(include=['object', 'category']).columns
        dummy_vars = pd.get_dummies(covariates[non_numeric_cols], drop_first=True)
        numeric_vars = covariates.drop(columns=non_numeric_cols)
        metadata_reformed = pd.concat([metadata[['SAMPLE']], numeric_vars, dummy_vars], axis=1)
        
        # Merge all data - merge genetic data, methylation data, and metadata on SAMPLE and HAPLOTYPE columns (genetics + methylation)
        # This is based on the assumption that for each sample, genetics H1 and methylation H1 match
        # In future offer command line option to examine all haplotype combinations per sample, or trans haplotype combos
        # merged data if phased
        if (args.simulate_unphased is False):
            merged_data = genetic_data.merge(methylation_data, on=['SAMPLE','HAPLOTYPE']).merge(metadata_reformed, on='SAMPLE')
        # merged data if unphased
        else:
            merged_data = unphased_genetic_data.merge(unphased_methylation_data, on=['SAMPLE']).merge(metadata_reformed, on='SAMPLE')
        
        # Subset to relevant columns
        # I think TARGET is the name of the methylation region
        outcomes = list(set(probes['TARGET']) & set(merged_data.columns))
        # And VARIANT is the name of the SV/SNV associated
        predictors = list(set(variants['VARIANT']) & set(merged_data.columns))
        covariates_list = list(set(metadata_reformed.columns) & set(merged_data.columns))
        # remove SAMPLE from covariates
        covariates_list.remove('SAMPLE')
        
        # Run regression models and collect results
        for outcome in outcomes:
            for predictor in predictors:
                try:
                    # different procedure depending on whether or not unphased data simulated
                    
                    if (args.simulate_unphased is True):
                        # reference frequency is proportion of 0 genotypes
                        ref_variant_frequency=sum(haplotype_outcome_predictor[predictor]==0)/len(haplotype_outcome_predictor[predictor])
                        # alternate frequency is proportion of 1 genotypes
                        alt_variant_frequency=sum(haplotype_outcome_predictor[predictor]==1)/len(haplotype_outcome_predictor[predictor])
                        # minor allele frequency is least of either reference or alternate frequencies
                        overall_minor_allele_frequency=min(ref_variant_frequency,alt_variant_frequency)
                        
                        # only run if simulated unphased minor allele frequency is over 0.05
                        if (overall_minor_allele_frequency > args.overall_MAF):
                            
                            formula = f"Q('{outcome}') ~ {predictor} + {' + '.join(covariates_list)}"
                            model = smf.ols(formula, data=merged_data).fit()
                            results.append({
                                'gene': gene_of_interest,
                                'window_size': args.window_size,
                                'outcome': outcome,
                                'predictor': predictor,
                                'beta': model.params.get(predictor, np.nan),
                                # add age coefficient
                                # add sex_male (binary) coefficient
                                # add PMI coefficient
                                'std_err': model.bse.get(predictor, np.nan),
                                'r2': model.rsquared,
                                'p_value': model.pvalues.get(predictor, np.nan),
                                'N': model.nobs,
                                'predictor_freq': merged_data[predictor].mean(skipna=True),
                                # print intercept
                                'intercept': model.params.Intercept,
                                # print overall allele frequency
                                'overall_MAF': overall_minor_allele_frequency,
                            })
                        
                    else:
                        # check for overall minor allele frequency (lesser of two absent or present frequencies)
                        haplotype_outcome_predictor=merged_data[['HAPLOTYPE',outcome,predictor]]
                        # reference frequency is proportion of 0 genotypes
                        ref_variant_frequency=sum(haplotype_outcome_predictor[predictor]==0)/len(haplotype_outcome_predictor[predictor])
                        # per haplotype
                        ref_variant_frequency_H1=sum(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][predictor]==0)/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][predictor])
                        ref_variant_frequency_H2=sum(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][predictor]==0)/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][predictor])
                        # alternate frequency is proportion of 1 genotypes
                        alt_variant_frequency=sum(haplotype_outcome_predictor[predictor]==1)/len(haplotype_outcome_predictor[predictor])
                        # per haplotype
                        alt_variant_frequency_H1=sum(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][predictor]==1)/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][predictor])
                        alt_variant_frequency_H2=sum(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][predictor]==1)/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][predictor])
                        # minor allele frequency is least of either reference or alternate frequencies
                        overall_minor_allele_frequency=min(ref_variant_frequency,alt_variant_frequency)
                        # MAF per haplotype
                        minor_allele_frequency_H1=min(ref_variant_frequency_H1,alt_variant_frequency_H1)
                        minor_allele_frequency_H2=min(ref_variant_frequency_H2,alt_variant_frequency_H2)
                        # get single haplotype minor allele frequency
                        # check for data with missing methylation for one haplotype
                        proportion_meth_H1_missing=haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][outcome].isna().sum()/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H1'][outcome])
                        proportion_meth_H2_missing=haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][outcome].isna().sum()/len(haplotype_outcome_predictor[haplotype_outcome_predictor['HAPLOTYPE']=='H2'][outcome])
                        
                        if (overall_minor_allele_frequency > args.overall_MAF and proportion_meth_H1_missing < args.per_haplotype_missing_methylation_rate and proportion_meth_H2_missing < args.per_haplotype_missing_methylation_rate) or (proportion_meth_H1_missing > args.per_haplotype_missing_methylation_rate and minor_allele_frequency_H2 > args.per_haplotype_MAF) or (proportion_meth_H2_missing > args.per_haplotype_missing_methylation_rate and minor_allele_frequency_H1 > args.per_haplotype_MAF):
                        
                            formula = f"Q('{outcome}') ~ {predictor} + {' + '.join(covariates_list)}"
                            model = smf.ols(formula, data=merged_data).fit()
                        
                            results.append({
                                'gene': gene_of_interest,
                                'window_size': args.window_size,
                                'outcome': outcome,
                                'predictor': predictor,
                                'beta': model.params.get(predictor, np.nan),
                                # add age coefficient
                                # add sex_male (binary) coefficient
                                # add PMI coefficient
                                'std_err': model.bse.get(predictor, np.nan),
                                'r2': model.rsquared,
                                'p_value': model.pvalues.get(predictor, np.nan),
                                'N': model.nobs,
                                'predictor_freq': merged_data[predictor].mean(skipna=True),
                                # print intercept
                                'intercept': model.params.Intercept,
                                # print overall allele frequency
                                'overall_MAF': overall_minor_allele_frequency,
                                # print minor allele frequency for haplotype 1 (H1)
                                'H1_MAF': minor_allele_frequency_H1,
                                # print minor allele frequency for haplotype 2 (H2)
                                'H2_MAF': minor_allele_frequency_H2,
                                # print missing methylation per haplotype
                                'H1_missing_meth': proportion_meth_H1_missing,
                                'H2_missing_meth': proportion_meth_H2_missing
                            })
                    
                    # print(model.summary())
                except Exception as e:
                    print(f"Regression failed for {outcome} ~ {predictor}: {e}")
    
    # Save results
    df_results = pd.DataFrame(results)
    df_results.to_csv(args.output, index=False)
    
    # Print runtime and max RAM usage
    end_time = time.time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    print(f"Max RAM Usage: {process.memory_info().rss / (1024 ** 2):.2f} MB")

if __name__ == "__main__":
    main()