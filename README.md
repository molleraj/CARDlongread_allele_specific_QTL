# NIA CARD Long Read Sequencing Allele-Specific Quantitative Trait Locus (QTL) Analysis Pipeline
This is the repository for the allele-specifc QTL pipeline developed by the CARD long-read sequencing/applied neurogenomics group. We developed a pipeline to perform allele-specific quantitative trait locus (QTL) analysis on phased, harmonized genetic variant and methylation data initially generated for NABEC and HBCC cohorts. Existing QTL pipelines take unphased variant and phenotype data as input. We anticipated that an allele-specific analysis might capture haplotype-specific signals averaged out in unphased data, particularly for methylation. Future development will extend this pipeline to handle more types of omic phenotypic data.

## Pipeline overview
## Dependencies
## Containerization
## Metadata considerations
## Variant filtering
We have initially filtered SNVs and SVs used for the allele-specific QTL analysis both to match filtering parameters used in [our past standard mQTL analyses](https://www.biorxiv.org/content/10.1101/2024.12.16.628723v1) and to reduce run time and significant variants to examine. Our default filtering parameters are variants with a minor allele frequency (amongst samples) of over 5%, a genotyping rate of 95% or higher, and a Hardy-Weinberg equilibrium p-value cutoff of 0.001. We performed additional, initial linkage disequilibrium (LD) pruning for our allele-specific mQTL, pruning variants with a pairwise correlation coefficient (r<sup>2</sup>) of 0.3 or higher within a 1000 variant count window with 50 variant steps between windows.
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
Input data for the QTL pipeline must be standardized in the formats listed in the [CARDlongread_data_standardization](https://github.com/NIH-CARD/CARDlongread_data_standardization) repository using the included scripts. The goal of this standardization is to convert metadata, phased genetic variants, and phased methylation data into a haplotype-specific form ready for downstream QTL and machine learning analysis. Filtered variants from the step above should be used as genetic variant input to generate respective maps and data matrices.

## Running the allele-specific QTL
```
usage: long_read_QTLs.py [-h] --chromosome CHROMOSOME --roi_map ROI_MAP --genetic_data GENETIC_DATA --methylation_data METHYLATION_DATA --genetic_map GENETIC_MAP --methylation_map METHYLATION_MAP --metadata METADATA --output OUTPUT
                         [--window_size WINDOW_SIZE]

Perform QTL analysis with user-specified parameters. Input data formats described in CARDlongread_data_standardization repository.

optional arguments:
  -h, --help            show this help message and exit
  --chromosome CHROMOSOME
                        Specify the chromosome to subset the ROI map (e.g., 'chr1').
  --roi_map ROI_MAP     Path to the ROI map file.
  --genetic_data GENETIC_DATA
                        Path to the genetic data file (e.g., '/path/to/genetic_data').
  --methylation_data METHYLATION_DATA
                        Path to the methylation data file.
  --genetic_map GENETIC_MAP
                        Path to the genetic map file.
  --methylation_map METHYLATION_MAP
                        Path to the methylation map file.
  --metadata METADATA   Path to the metadata file.
  --output OUTPUT       Output file destination for results.
  --window_size WINDOW_SIZE
                        Window size for defining gene regions (default: 500000).
```
## Postprocessing and data visualization
```
usage: long_read_QTL_lambda.py [-h] --input INPUT

Perform lambda (genomic inflation factor) calculation for QTL results.

optional arguments:
  -h, --help     show this help message and exit
  --input INPUT  Path to input QTL file for all chromosomes; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.
```
```
usage: long_read_QTL_fdr_correction.py [-h] --input INPUT --output OUTPUT [--rejected | --no-rejected] [--drop_na | --no-drop_na] [--top_hit_per_gene | --no-top_hit_per_gene] [--top_hit_per_region | --no-top_hit_per_region]
                                       [--unique_hits | --no-unique_hits]

Perform multiple hypothesis testing false-discovery rate (FDR) correction for allele-specific methylation QTL output from all chromosomes.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Path to input QTL file for all chromosomes; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.
  --output OUTPUT       Path to the output QTL file with additional two fields - rejected and pvalue-corrected.
  --rejected, --no-rejected
                        Set this option to only print results where the null hypothesis was rejected (i.e., statistically significant associations). (default: False)
  --drop_na, --no-drop_na
                        Drop any cases where p value is listed as NA. (default: False)
  --top_hit_per_gene, --no-top_hit_per_gene
                        Print only most significant hit per gene. (default: False)
  --top_hit_per_region, --no-top_hit_per_region
                        Print only most significant hit per methylation region. (default: False)
  --unique_hits, --no-unique_hits
                        Print only unique set of included tests. (default: False)
```
```
usage: long_read_QTL_data_visualization.py [-h] --input_combinations INPUT_COMBINATIONS --genetic_data GENETIC_DATA --methylation_data METHYLATION_DATA --output_prefix OUTPUT_PREFIX [--trans_comparison | --no-trans_comparison]
                                           [--violin_plot | --no-violin_plot] [--strip_plot | --no-strip_plot]

Generate boxplot/violinplot phenotype distribution visualizations and merged tables for particular phenotype/variant combinations.

optional arguments:
  -h, --help            show this help message and exit
  --input_combinations INPUT_COMBINATIONS
                        Headerless CSV file with list of comma-separated phenotype/variant combinations to be analyzed.
  --genetic_data GENETIC_DATA
                        Path to genetic data input file; format described in CARDlongread_data_standardization repository.
  --methylation_data METHYLATION_DATA
                        Path to methylation data input file; format described in CARDlongread_data_standardization repository.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output files (plot and table per variant/phenotype combination).
  --trans_comparison, --no-trans_comparison
                        Show violin plots instead of box plots (optional; default false) (default: False)
  --violin_plot, --no-violin_plot
                        Show violin plots instead of box plots (optional; default false) (default: False)
  --strip_plot, --no-strip_plot
                        Show strip plots instead of swarm plots inside box/violin plots (optional; default false) (default: False)
```
## Comparison with standard QTL results
```
usage: long_read_QTL_comparison.py [-h] --tensor_QTL TENSOR_QTL --allele_specific_QTL ALLELE_SPECIFIC_QTL --region_type REGION_TYPE --output_prefix OUTPUT_PREFIX

Compare allele-specific QTL output to tensorQTL output on common methylation region. Print common regions and those unique to each.

options:
  -h, --help            show this help message and exit
  --tensor_QTL TENSOR_QTL
                        Path to input tensorQTL tsv file for comparisons; 'phenotype_id variant_id af pval_nominal slope slope_se pval_perm bh_fdr qval Chromosome TOP SV ID TOP SV Causal Post Probablity TOP SNV ID TOP SNV Causal Post
                        Probablity' as header.
  --allele_specific_QTL ALLELE_SPECIFIC_QTL
                        Path to input allele-specific QTL file for comparisons; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.
  --region_type REGION_TYPE
                        Region type (CGI for CpG islands, GB for gene bodies, PROM for promoters).
  --output_prefix OUTPUT_PREFIX
                        Prefix for output files with tensorQTL only, allele-specific QTL only, and common hits.
```
