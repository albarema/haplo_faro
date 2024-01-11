# Chromosome painting Faroese

Population structure and local ancestry inference of a certain population (we used it for Faroese's) using a window-based (phased) haplotype clustering approach. The neural network (NN) can be used in whole-genome data. Here are more details about the neural network (NN): [Haplonet github](https://github.com/Rosemeis/HaploNet).

Software requirements:
- Snakemake 
- Conda
- bcftools 
- Haplonet (which relies on other Python packages, follow instructions [here](https://github.com/Rosemeis/HaploNet?tab=readme-ov-file#install-and-build))
- R (version >= 4.0)
- Python (>= 3.6)

Please read [this guide](guide_snakemake_conda_env.md) if you have questions about Snakemake and conda environments (corresponding to the file ```guide_snakemake_conda_env.md```). 

## Getting started

1. Clone this repository (git clone XXXX)
2. Create a conda environment using the environmental.yaml file provided (conda env create --file environment.yaml)
3. Prepare all input files

Input files:
- phased (and imputed) whole-genome sequencing (WGS) data in VCF/BCF format (NB: ideally one file per chromosome, otherwise, check what to do when using a merged file ```--like```)
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
First of all, we will run this file [step1_training_popst](rules/haplonet_main.smk). 

```
snakemake --snakefile rules/haplonet_main.smk -j10
```
Information on what each rule does in the snakemake file but basically, it first outputs the NN log-likelihoods which are later used to infer global population structure (PCA and admixture). We will use 10 seeds (more seeds are required for higher Ks). For the chosen K, check if the log-likelihoods have converged, otherwise, run up to 100 seeds. 

The way we estimate the admixture proportions is by using a fixed "F" matrix estimated with other samples which are considered as a "training set" (important not to include the population of interest, Faroese in this case). Then we estimate the "Q" matrix for all samples. 

Output files:
- log-likelihoods in binary NumPy matrix per chromosome (*.loglike.npy)
- log file with parameters used in the training (*.log)
- ancestry proportions in a text file (*.q)
- ancestral cluster frequencies in a binary NumPy matrix (*.f.npy)
- eigenvectors directly using singular value decomposition (SVD) (*.eigenvecs)
  
Then, use the admixture results with the highest log-likelihood for further analyses. The Snakemake file [step2_best_plotting](rules/plotting_haplonet.smk) does this for you. 
 
```
snakemake --snakefile rules/plotting_haplonet.smk -j10
```

Finally, we need to run [step3_fastash](rules/fastash.smk) to get chromosome painting:
```
snakemake --snakefile rules/fatash.smk -j10
```
Output file:
- Best cluster per window (*.path, window-size set in step 1). Each line corresponds to a haplotype. Haplotypes 1 and 2 from the same individual are found consecutive and the order of the individuals is the same as in the VCF file. 

