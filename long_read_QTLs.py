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
    parser.add_argument("--window_size", type=int, default=500000, help="Window size for defining gene regions (default: 500000).")
    return parser.parse_args()

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
        merged_data = genetic_data.merge(methylation_data, on=['SAMPLE','HAPLOTYPE']).merge(metadata_reformed, on='SAMPLE')
        
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
                    formula = f"Q('{outcome}') ~ {predictor} + {' + '.join(covariates_list)}"
                    model = smf.ols(formula, data=merged_data).fit()
                    
                    results.append({
                        'gene': gene_of_interest,
                        'window_size': args.window_size,
                        'outcome': outcome,
                        'predictor': predictor,
                        'beta': model.params.get(predictor, np.nan),
                        'std_err': model.bse.get(predictor, np.nan),
                        'r2': model.rsquared,
                        'p_value': model.pvalues.get(predictor, np.nan),
                        'N': model.nobs,
                        'predictor_freq': merged_data[predictor].mean(skipna=True)
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
