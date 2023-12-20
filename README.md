# haplo_faro

Population structure and local ancestry inference of a certain population (we used it for Faroese's) using a window-based (phased) haplotype clustering approach. The neural network (NN) can be used in whole-genome data.  Here are more details about the NN: [Haplonet github](https://github.com/Rosemeis/HaploNet).

Software requirements:
- Snakemake 
- Conda
- bcftools 
- Haplonet
- R (version >= 4.0)
- Python (>= 3.6)

Please read the ```guide_snakemake_conda_env.md``` if you have questions about snakemake and conda environments. 

Input:
- phased (and imputed) data (VCF format)
- tab-delimited file with the sample IDs (same order as VCF)

Pre-step:  since we aimed to paint the Faroese chromosomes using ancient samples we merged both datasets using ```rules/prep_merged_vcf.smk```. You can skip this step if all samples are already in the same VCF file. 
First step: filter the VCF using bcftools to remove missing data and to only keep bi-allelic sites with a MAF > 0.05.
```
bcftools view -m 2 -M 2 -v snps -i 'MAF > 0.05'
```

Run Haplonet to infer population structure and local ancestry using the Snakemake file provided. 

Since 10 seeds have been used, use the admixture results with the highest log-likelihood for further analyses. The Snakemake file ```plotting_haplonet.smk``` does this for you.  
