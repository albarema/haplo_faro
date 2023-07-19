# haplo_faro

Population structure and local ancestry inference of a certain population (we used it for Faroese's) using a window-based (phased) haplotype clustering approach. The neural network (NN) can be used in whole-genome data.  More details about the NN here: [Haplonet github](https://github.com/Rosemeis/HaploNet).

Requierements: phased data. 

First step: filter the VCF using bcftools to remove missing data and to only keep bi-allelic sites with a MAF > 0.05. 
```
bcftools view -m 2 -M 2 -v snps -i 'MAF > 0.05'
```

Run haplonet to infer population structure and local ancestry using the Snakemake file provided. 

Since 10 seeds have been used, use the admixture results with the highest log-likelihood for further analyses. The Snakemake file plotting_haplonet.smk does this for you.  
