# NIA CARD Long Read Sequencing Allele-Specific Quantitative Trait Locus (QTL) Analysis Pipeline
This is the repository for the allele-specifc QTL pipeline developed by the CARD long-read sequencing/applied neurogenomics group. We developed a pipeline to perform allele-specific quantitative trait locus (QTL) analysis on phased, harmonized genetic variant and methylation data initially generated for NABEC and HBCC cohorts. Existing QTL pipelines take unphased variant and phenotype data as input. We anticipated that an allele-specific analysis might capture haplotype-specific signals washed out in unphased data, particularly for methylation. Future development will extend this pipeline to handle more types of omic phenotypic data.
## Pipeline overview
## Dependencies
## Containerization
## Variant filtering
We have initially filtered SNVs and SVs used for the allele-specific QTL analysis both to match filtering parameters used in [our past standard mQTL analyses](https://www.biorxiv.org/content/10.1101/2024.12.16.628723v1) and to reduce run time and significant variants to examine.
```
Usage: ./long_read_QTL_initial_variant_filtering.sh -i input_prefix -s sample_exclude_list -m maf_cutoff -g missing_genotype_rate -h hwe_pvalue -p indep_pairwise_ld_pruning_values
	-i Input VCF file prefix
	-s Path to list of samples to exclude
	-m Minor allele frequency (MAF) cutoff (default 0.05)
	-g Missing genotype rate cutoff (default 0.05)
	-h Hardy-Weinberg equilibrium p-value cutoff (default 0.001)
	-p Independent pairwise LD pruning settings (default "1000 50 0.3"). Numbers indicate window, step size, and r2 value for indep-pairwise pruning.
```
## Input data standardization
Input data for the QTL pipeline must be standardized in the formats listed in the [CARDlongread_data_standardization](https://github.com/NIH-CARD/CARDlongread_data_standardization) repository using the included scripts. 
## Running the allele-specific QTL
## Postprocessing and data visualization
## Comparison with standard QTL results
