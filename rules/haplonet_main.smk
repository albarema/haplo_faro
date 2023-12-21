# run from /home/alba/haplo
configfile: 'config_haplo.yml'
import random 

#-----------------------------------------------------------------------------------------------
CHR=range(1,23)
#SEEDN = random.sample(range(0, 1000), 10)
SEEDN=range(1,101)
KS=range(3,7)
KS_super = range(5,9)

wildcard_constraints:
    chrn= "[0-9]+"

#-----------------------------------------------------------------------------------------------
rule all:
    input:
        # Training NN to get haplotype clusters
        # expand("{panel}/chrs-{windowSize}.loglike.list", panel=['ancient.3354BP'],windowSize='ws'),
        # admixture Q estimates
        # expand("{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.q", panel=['ancient.3354BP'],windowSize='ws',ks=6, seedN=SEEDN),
        # concatenate all admixture results and get the best seed run
        # expand("{panel}/allbest-unsuperv-{windowSize}.ks-admixture.txt", panel=['ancient.3354BP'],windowSize='wl'), 
        # unsupervised-supervised haplonet 
        # expand("{panel}-unsupervk5/haplonet-superv-{windowSize}.admixture.k{ks}.s{seedN}.q", panel=['ancient.3354BP'],windowSize='ws',ks=KS_super, seedN=SEEDN),


rule training:
    "training NN haplotype clustering"
    input:
        os.path.join(config['in_path'],"vcf_chrs_ancient/FOFG-WGS-MESO3354BP-chr{chrn}.vcf.gz"),
    output:
        out="{panel}/chr{chrn}-{windowSize}.loglike.npy", 
        args=temp("{panel}/chr{chrn}-{windowSize}.args")
    log: "{panel}/chr{chrn}-{windowSize}.log"
    params: 
        wsize=lambda wildcards: 512 if wildcards.windowSize == "ws" else (1024 if wildcards.windowSize == "wl" else 1024),
        out="{panel}/chr{chrn}-{windowSize}"
    threads: 5
    shell:
        "haplonet train --vcf {input} "
        "-t {threads} "
        "--subsplit 2 "
        "--x_dim {params.wsize} "
        "--out {params.out} 2> {log}"

rule filepath:
    "concat all chroms"
    input:
        expand("{panel}/chr{chrn}-{windowSize}.loglike.npy", chrn=CHR, allow_missing=True)
    output:
        "{panel}/chrs-{windowSize}.loglike.list"
    shell:
        "realpath {input} > {output}"

rule get_clusters:
    "supervised admixture proportion assignments and metadata completion"
    input:
        vorder="samples/MESO3354BP_vcf_order.txt",
        meta="/home/alba/mesoneo/neo.impute.1000g.sampleInfo_clusterInfo.txt",
    output:
        super="samples/MESO3354BP.supervised.clusters.txt",
        unsuper="samples/MESO3354BP.training.clusters.txt",
        meta="metadata/MESO3354BP.sampleInfo_clusterInfo.txt"
    shell:
        "Rscript cluster_ks_haplo.r "
        "-m {input.meta} "
        "-s {input.vorder} "
        "-o {output.meta} "
        "-t {output.unsuper} "
        "-v {output.super}"


rule unsuper:
    "admixture unsupervised-version, update q for Faroese data"
    input:
        loglist="{panel}/chrs-{windowSize}.loglike.list",
        training="samples/MESO3354BP.training.clusters.txt",
    output:
        out="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.q"
    params:
      out="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}"
    log: "{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.log"
    threads: 30
    shell:
        "haplonet admix --filelist {input.loglist} "
        "--K {wildcards.ks} "
        "--threads {threads} "
        "--training {input.training} "
        "--seed {wildcards.seedN} "
        "--out {params.out} > {log}"

rule best_ks:
    params: "{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}"
    output:
        tmp=temp("{panel}/best-unsuperv-{windowSize}.k{ks}.txt")
    shell:
        """
        cat <(paste -d"\t" <(head -n1 --quiet {params}*.log | awk '{{print $8}}') \
        <(grep -m 1 "Final" {params}*.log | nl | awk -v OFS='\t' '{{print {wildcards.ks},$4}}'))| \
        sort -k 3 | head -n1 > {output.tmp}
        """
  
rule allBest:
    input:
        expand("{panel}/best-unsuperv-{windowSize}.k{ks}.txt", ks=KS, allow_missing=True)
    output:
        "{panel}/allbest-unsuperv-{windowSize}.ks-admixture.txt"
    shell:
        """
        cat <(echo -e "seed\tk\tlog-likelihood") <(cat {input} | sed -e 's/seed=//' ) > {output}     
        """

rule super:
    "admixture supervised-version, using as sources individuals with >.99 from the unsupervised-version"
    input:
        loglist="{panel}/chrs-{windowSize}.loglike.list",
        supervised="samples/ANC3354BP.unsupervk5.supervised.clusters.txt",
    output:
        out="{panel}-unsupervk5/haplonet-superv-{windowSize}.admixture.k{ks}.s{seedN}.q"
    params:
        out="{panel}-unsupervk5/haplonet-superv-{windowSize}.admixture.k{ks}.s{seedN}"
    log: "{panel}-unsupervk5/haplonet-superv-{windowSize}.admixture.k{ks}.s{seedN}.log"
    threads: 30
    shell:
        "haplonet admix --filelist {input.loglist} "
        "--K {wildcards.ks} "
        "--threads {threads} "
        "--supervised {input.supervised} "
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
