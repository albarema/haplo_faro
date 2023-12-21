# Failed to solve scheduling problem with ILP solver. Falling back to greedy solver. Run Snakemake with --verbose to see the full solver output for debugging the problem. snakemake --snakefile ../rules/haplonet.smk -n
# run from /home/alba/haplo
configfile: 'config_haplo.yml'
import random 

#-----------------------------------------------------------------------------------------------
CHR=range(1,23)
SEEDN = random.sample(range(0, 1000), 10)
KS=range(2,9)

#-----------------------------------------------------------------------------------------------

rule all:
    input:
        expand("{panel}/chrs-{windowSize}.loglike.list", panel=['ancient'],windowSize='wl'),
        expand("{panel}/haplonet-{windowSize}.admixture.k{ks}.s{seedN}.q", panel=['ancient'],windowSize='wl',ks=KS, seedN=SEEDN),
        expand("{panel}/haplonet-{windowSize}.eigenvecs", panel=['ancient'],windowSize='wl')
        

rule training:
    input:
        os.path.join(config['in_path'],"vcf_chrs_{panel}/FOFG-WGS-MESO-chr{chrn}.vcf.gz"),
    output:
        out="{panel}/chr{chrn}-{windowSize}.loglike.npy", 
        args=temp("{panel}/chr{chrn}-{windowSize}.args")
    log: "{panel}/chr{chrn}-{windowSize}.log"
    params: 
        wsize=lambda wildcards: 512 if wildcards.windowSize == "ws" else (1024 if wildcards.windowSize == "wl" else 1024),
        out="{panel}/chr{chrn}-{windowSize}"
    threads: 15
    shell:
        "haplonet train --vcf {input} "
        "-t {threads} "
        "--x_dim {params.wsize} "
        "--out {params.out} 2> {log}"

rule filepath:
    input:
        expand("{panel}/chr{chrn}-{windowSize}.loglike.npy", chrn=CHR, allow_missing=True)
    output:
        "{panel}/chrs-{windowSize}.loglike.list"
    shell:
        "realpath {input} > {output}"


rule ancestry:
    input:
        "{panel}/chrs-{windowSize}.loglike.list"
    output:
        "{panel}/haplonet-{windowSize}.admixture.k{ks}.s{seedN}.q"
    params:
        out="{panel}/haplonet-{windowSize}.admixture.k{ks}.s{seedN}"
    log: "{panel}/{windowSize}.admixture.k{ks}.s{seedN}.log"
    threads: 30
    shell:
        "haplonet admix --filelist {input} "
        "--K {wildcards.ks} "
        "--threads {threads} "
        "--seed {wildcards.seedN} "
        "--out {params.out} > {log}"

rule pca:
    input:
        "{panel}/chrs-{windowSize}.loglike.list"
    output:
        out="{panel}/haplonet-{windowSize}.eigenvecs",
    params:
        out="{panel}/haplonet-{windowSize}"
    log: "{panel}/pca.haplonet-{windowSize}.log"
    threads: 15
    shell:
        "haplonet pca --filelist {input} "
        "--threads {threads} "
        "--out {params.out} > {log}"
  
