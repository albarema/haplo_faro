# Chromosome painting Faroese

Population structure and local ancestry inference of a certain population (we used it for Faroese's) using a window-based (phased) haplotype clustering approach. The neural network (NN) can be used in whole-genome data. Here are more details about the NN: [Haplonet github](https://github.com/Rosemeis/HaploNet).

Software requirements:
- Snakemake 
- Conda
- bcftools 
- Haplonet
- R (version >= 4.0)
- Python (>= 3.6)

Please read [this guide](guide_snakemake_conda_env.md) if you have questions about Snakemake and conda environments (```guide_snakemake_conda_env.md```). 

Input files:
- phased (and imputed) data (VCF format)
- tab-delimited file with the sample IDs (same order as VCF)

## Workflow overview
- Pre-step: combine datasets (optional)
- Pre-processing: variant- and individual-level
- Haplonet
    - Training (NN)
    - Admixture (training + supervised)
    - Fatash
    - PCA
- Plotting

1. Pre-step:  since we aimed to paint the Faroese chromosomes using ancient samples we merged both datasets using ```rules/prep_merged_vcf.smk```. You can skip this step if all samples are already in the same VCF file. 
2. Pre-processing: filter the VCF using bcftools to remove missing data and to only keep bi-allelic sites with MAF > 0.05.
```
bcftools view -m 2 -M 2 -v snps -i 'MAF > 0.05'
```
3. Run Haplonet to infer population structure and local ancestry using the Snakemake file provided. 

We will be running 10 seeds (more seeds are required for higher Ks). For the chosen K, check if the log-likelihoods have converged, otherwise, run up to 100 seeds. Then, use the admixture results with the highest log-likelihood for further analyses. The Snakemake file ```plotting_haplonet.smk``` does this for you.  
