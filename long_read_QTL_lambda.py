# calculate lambda value (genomic inflation factor) for full QTL results

import pandas as pd
import numpy as np
import argparse
import os
from scipy.stats import chi2

def parse_args():
    parser = argparse.ArgumentParser(description="Perform lambda (genomic inflation factor) calculation for QTL results.")
    parser.add_argument("--input", required=True, help="Path to input QTL file for all chromosomes; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.")
    return parser.parse_args()

def main():
    # Parse input and output arguments.
    args = parse_args()
    # load input QTL table
    input_QTL_df = pd.read_csv(args.input)
    # must drop na to run...
    input_QTL_df = input_QTL_df[input_QTL_df.p_value.isna() != True]

    # Get p values from QTL table
    p_values = input_QTL_df['p_value']

    # Convert p-values to chi-square statistics with 1 degree of freedom
    chi2_stats = chi2.isf(p_values, df=1)  # inverse survival function

    # Lambda = median of observed chi2 / median of expected chi2 (df=1)
    median_chi2_obs = np.median(chi2_stats)
    median_chi2_null = chi2.ppf(0.5, df=1)  # ~0.455
    lambda_gc = median_chi2_obs / median_chi2_null

    # Print results
    print(f"Genomic Inflation Factor (λ): {lambda_gc:.4f}")

if __name__ == "__main__":
    main()