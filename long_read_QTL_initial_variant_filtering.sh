#!/bin/bash
# initial script for cleaning up input variants, whether SV or SNV
# load necessary modules
module load plink bcftools
# print help statements
helpFunction()
{
   echo "Script for initial variant cleanup pre-allele specific QTL, using plink and bcftools. Splits multiallelic variants by default."
   echo "Usage: $0 -i input_prefix -s sample_exclude_list -m maf_cutoff -g missing_genotype_rate -h hwe_pvalue -p indep_pairwise_ld_pruning_values"
   echo -e "\t-i Input VCF file prefix"
   echo -e "\t-s Path to list of samples to exclude"
   echo -e "\t-m Minor allele frequency (MAF) cutoff (default 0.05)"
   echo -e "\t-g Missing genotype rate cutoff (default 0.05)"
   echo -e "\t-h Hardy-Weinberg equilibrium p-value cutoff (default 0.001)"
   echo -e "\t-p Independent pairwise LD pruning settings (default \"1000 50 0.3\"). Numbers indicate window, step size, and r2 value for indep-pairwise pruning."
   exit 1 # Exit script after printing help
}

# set default values
maf_cutoff=0.05
missing_genotype_rate=0.05
hwe_pvalue=0.001
indep_pairwise_ld_pruning_values="1000 50 0.3"

# get options listed above
# defaults overridden if necessary
while getopts "i:s:m:g:h:p:" opt
do
   case "$opt" in
      i ) input_prefix="$OPTARG" ;;
      s ) sample_exclude_list="$OPTARG" ;;
      m ) maf_cutoff="$OPTARG" ;;
      g ) missing_genotype_rate="$OPTARG" ;;
      h ) hwe_pvalue="$OPTARG" ;;
      p ) indep_pairwise_ld_pruning_values="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case required parameters are empty
if [ -z "$input_prefix" ]
then
   echo "The required parameter (-i input_prefix) is empty.";
   helpFunction
fi

# different chained commands depending on whether excluded sample list provided or not
# in both cases, split all multiallelic variants, filter based on MAF/genotype/HWE, and then LD prune with indep-pairwise argument
if [ -z "$sample_exclude_list" ]
then
    	bcftools norm -m -both $input_prefix.vcf.gz -O z -o $input_prefix"_"multsplit.vcf.gz
        plink2 --vcf $input_prefix"_"multsplit.vcf.gz --output-chr chrMT --maf $maf_cutoff --geno $missing_genotype_rate --hwe $hwe_pvalue --vcf-half-call m --rm-dup force-first --export vcf bgz --out $input_prefix"_"multsplit_filtered
        plink2 --vcf $input_prefix"_"multsplit_filtered.vcf.gz --output-chr chrMT --indep-pairwise $indep_pairwise_ld_pruning_values --out $input_prefix"_"multsplit_filtered_pruned
	# extract variants listed as pruned 
	plink2 --extract $input_prefix"_"multsplit_filtered_pruned.prune.in --output-chr chrMT --vcf $input_prefix"_"multsplit_filtered.vcf.gz --export vcf bgz --out $input_prefix"_"multsplit_filtered_pruned
else
	# first exclude samples with bcftools view -S ^exclude_list_filename
	bcftools view -S "^"$sample_exclude_list $input_prefix.vcf.gz -O z -o $input_prefix"_"sample_excluded.vcf.gz
	bcftools norm -m -both $input_prefix"_"sample_excluded.vcf.gz -O z -o $input_prefix"_"sample_excluded_multsplit.vcf.gz
	plink2 --vcf $input_prefix"_"sample_excluded_multsplit.vcf.gz --output-chr chrMT --maf $maf_cutoff --geno $missing_genotype_rate --hwe $hwe_pvalue --vcf-half-call m --rm-dup force-first --export vcf bgz --out $input_prefix"_"sample_excluded_multsplit_filtered
	plink2 --vcf $input_prefix"_"sample_excluded_multsplit_filtered.vcf.gz --output-chr chrMT --indep-pairwise $indep_pairwise_ld_pruning_values --out $input_prefix"_"sample_excluded_multsplit_filtered_pruned
	plink2 --extract $input_prefix"_"sample_excluded_multsplit_filtered_pruned.prune.in --vcf $input_prefix"_"sample_excluded_multsplit_filtered.vcf.gz --output-chr chrMT --export vcf bgz --out $input_prefix"_"sample_excluded_multsplit_filtered_pruned
fi
