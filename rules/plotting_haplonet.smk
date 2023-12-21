# snakemake --snakefile ../rules/plotting_haplonet.smk -n
# run from /home/alba/haplo
configfile: 'config_haplo.yml'

import pandas as pd


#-----------------------------------------------------------------------------------------------
KS=range(2,9)
CHR=range(1,23)


def filepath_best(wildcards):
    d = pd.read_csv((wildcards.panel)+"/allbest.ks-admixture.txt", sep="\t")
    print("d")
    seedv = str(d.loc[d['k'] == int(wildcards.ks)]['seed'][0])
    filepathB = (wildcards.panel)+"/haplonet.admixture.k"+ str(wildcards.ks)+ ".s" + seedv +".q"
    return filepathB

def nonbest_rm(wildcards):
    d = pd.read_csv(wildcards.panel+"/allbest-"+(wildcards.windowSize) + ".ks-admixture.txt", sep="\t")
    KS = str(wildcards.ks)
    seedv = str(d.loc[d['k'] == int(KS)]['seed'].values[0])
    filepathB = wildcards.panel+ "/*"+ wildcards.windowSize + ".admixture.k"+ str(wildcards.ks)+ ".s" + str(seedv) + "*"
    return filepathB
#-----------------------------------------------------------------------------------------------

rule all:
    input:
        # only one k
        # expand("{panel}/logs-{training}-{windowSize}.k{ks}.txt",training=['unsuperv'], panel=['ancient.3354BP'],windowSize=["ws","wl"], ks=KS)
        # needs to be run first
        # expand("{panel}/allbest-{training}-{windowSize}.ks-admixture.txt", training=['unsuperv'], panel=['ancient.3354BP'],windowSize="ws"), #'ancient.AtSc'
        # then: 
        # expand("{panel}/files2rm-{model}-{windowSize}_k{ks}.txt", panel=['ancient'], ks =KS,windowSize="ws")
        # expand("{panel}/test.{ks}.txt", panel=['ancient'], ks=2)
 
rule best_ks:
    params: "{panel}/haplonet-{training}-{windowSize}.admixture.k{ks}"
    output:
        tmp="{panel}/logs-{training}-{windowSize}.k{ks}.txt"
    shell:
        """
        cat <(paste -d"\t" <(head -n1 --quiet {params}*.log | awk '{{print $8}}') \
        <(grep -m 1 "Final" {params}*.log | nl | awk -v OFS='\t' '{{print {wildcards.ks},$4}}'))| \
        sort -k 3 > {output.tmp}
        """
  
rule allBest:
    input:
        expand("{panel}/logs-{training}-{windowSize}.k{ks}.txt", ks=KS, allow_missing=True)
    output:
        "{panel}/allbest-{training}-{windowSize}.ks-admixture.txt"
    shell:
        """
        cat <(echo -e "seed\tk\tlog-likelihood") <(head -n1 {input} | sed -e 's/seed=//' ) > {output}     
        """

rule best_ls:
    input: 
        a="{panel}/best-{model}-{windowSize}.k{ks}.txt",
        b="{panel}/allbest-{model}-{windowSize}.ks-admixture.txt"
    params: 
        fp= lambda wildcards: nonbest_rm(wildcards),
        lsp = "{panel}/*-{model}-{windowSize}.admixture.k{ks}*"
    output: 
        rmfile = "{panel}/files2rm-{model}-{windowSize}_k{ks}.txt"
    shell:
        """
        comm -2 -3 <(ls {params.lsp} ) <(ls {params.fp}| sort) | tail +2 > {output.rmfile} ;
        # rm $(cat output.rmfile)
        """

# Run this and use the sample for the plotting
# bcftools query -l /data/vcf_chrs_ancient/FOFG-WGS-MESO3600BP-chr1.vcf.gz > /home/alba/samples/ancient.3600BP_samples_vcf_order.txt

rule plota:
    input: filepath_best
    output: 
        "{panel}/figures/admixture.k{ks}.png"
    shell: 
        "Rscript scripts/plotting.haplonet.r {input}"

        
