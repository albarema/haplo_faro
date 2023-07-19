# snakemake --snakefile ../rules/plotting_haplonet.smk -n
# run from /home/alba/haplo
configfile: 'config_haplo.yml'

import pandas as pd


#-----------------------------------------------------------------------------------------------
KS=range(2,9)
CHR=range(1,22)


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
        # needs to be run first
        #expand("{panel}/allbest-{windowSize}.ks-admixture.txt", panel=['ancient'],windowSize="ws"),
        # then: 
        expand("{panel}/files2rm-{windowSize}_k{ks}.txt", panel=['ancient'], ks =KS,windowSize="ws")
        #expand("{panel}/test.{ks}.txt", panel=['ancient'], ks=2)
     
 
rule best_ks:
    params: "{panel}/{windowSize}.admixture.k{ks}"
    output:
        tmp=temp("{panel}/best-{windowSize}.k{ks}.tmp.txt"),
        tmp2=temp("{panel}/best-{windowSize}.k{ks}.txt")
    shell:
        """
        paste -d"\t" <(head -n1 --quiet {params}*.log | awk '{{print $8}}') \
        <(tail {params}*.log | grep "Final" | nl | awk -v OFS='\t' '{{print {wildcards.ks},$4}}') > {output.tmp};
        cat {output.tmp} | sort -r -k 3 | head -n1 > {output.tmp2}
        """
      
rule allBest:
    input:
        expand("{panel}/best-{windowSize}.k{ks}.txt", ks=KS, allow_missing=True)
    output:
        "{panel}/allbest-{windowSize}.ks-admixture.txt"
    shell:
        """
        cat <(echo -e "seed\tk\tlog-likelihood") <(cat {input} | sed -e 's/seed=//' ) > {output}     
        """

rule best_ls:
    input: 
        a="{panel}/best-{windowSize}.k{ks}.txt",
        b="{panel}/allbest-{windowSize}.ks-admixture.txt"
    params: 
        fp= lambda wildcards: nonbest_rm(wildcards),
        lsp = "{panel}/*{windowSize}.admixture.k{ks}*"
    output: 
        rmfile = "{panel}/torm/files2rm-{windowSize}_k{ks}.txt"
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
        