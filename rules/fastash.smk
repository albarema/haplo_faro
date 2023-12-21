# snakemake --snakefile ..rules/haplonet_fatash.smk
configfile: 'config_haplo.yml'

import pandas as pd


#-----------------------------------------------------------------------------------------------
KS=range(2,9)
CHR=range(1,23)

def get_vars(wildcards):
    d = pd.read_csv("samples/MESO3354BP_vcf_order.txt",names=["ids"])
    index_value_tuples = [(index*2-1, row['ids']) for index, row in d.iterrows()]
    index_int = d.index[d['ids'] == wildcards.sampleid].item()
    return index_value_tuples, index_int*2
#-----------------------------------------------------------------------------------------------


rule all:
    input:
        # estimate local ancestry
        # expand("{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.path", panel=['ancient.3354BP'],windowSize='ws',ks=5, seedN=91),
        # get painting files per sample 
        expand("{panel}/painting/allinds-unsuperv-{windowSize}.k{ks}.s{seedN}.tsv",panel=['ancient.3354BP'],windowSize='ws',ks=5, seedN=91)
     
 
rule local:
    input:
        loglist = "{panel}/chrs-{windowSize}.loglike.list",
        q="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.q",
        f="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.f.npy"
    output:
        prob="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.chr1.prob.npy",
        path="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.chr1.path",
    params:
        out="{panel}/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}"
    log: "{panel}/unsuperv-{windowSize}.fatash.k{ks}.s{seedN}.log"
    threads: 8
    wildcard_constraints: windowSize='[^-]*'
    shell:
        "haplonet fatash --filelist {input.loglist} "
        "--prop {input.q} "
        "--freq {input.f} "
        "--optim "
        "--threads {threads} "
        "--out {params.out} > {log}"

rule fatash_ind:
    input:
        ids="samples/MESO3354BP_vcf_order.txt",
        fatask=expand("ancient.3354BP/haplonet-unsuperv-{windowSize}.admixture.k{ks}.s{seedN}.chr{chrs}.path", chrs = CHR,allow_missing=True)
    output:
        painted="{panel}/painting/allinds-unsuperv-{windowSize}.k{ks}.s{seedN}.tsv"
    run:
        import pandas as pd
        
        # Create an empty DataFrame to store the combined data
        combined_data = pd.DataFrame()
        
        for chrs in range(1, 23):
            d = pd.read_csv("samples/MESO3354BP_vcf_order.txt",names=["ids"])
            input_file=pd.read_csv("ancient.3354BP/haplonet-unsuperv-ws.admixture.k5.s91.chr{}.path".format(chrs), names=["painting"]) 
        
            # Replicating each ID and maintaining the order
            replicated_ids = d['ids'].repeat(2).reset_index(drop=True)
            # creating haplotype order
            haplo_sequence = ['h1', 'h2'] * (len(replicated_ids) // 2)
        
            # creatind dataframe
            combined_df = pd.concat([replicated_ids, input_file], axis=1)
            combined_df['haplo'] = haplo_sequence
            combined_df['chromosome'] = chrs
        
            # Append the current DataFrame to the overall combined_data
            combined_data = pd.concat([combined_data, combined_df])
        
        combined_data.to_csv(output.painted, sep='\t', index=False) 
        