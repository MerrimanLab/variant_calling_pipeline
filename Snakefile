import pandas as pd
configfile: "config.yaml"


BUILDS= ["b37","b38"]


    
CHRS = list(range(1,23))

#b37_exome = config["b37_exome"]
#b38_exome = config["b38_exome"]

wildcard_constraints:
    build = "b3[78]"


nice_cmd="nice -n 10"


sequences = pd.read_table(config["sequences_metadata"], dtype=str).set_index(["unique_id"], drop=False)
sequences["run_info"] = sequences.unique_id.str.replace(r'^([A-Za-z\d-]+_)', "")


samplesInGroups = {}
for sample in sequences["Sample"].unique():
  samplesInGroups[sample] = sequences.loc[(sequences['Sample'] == sample), ["unique_id"]]["unique_id"].to_list()

def get_build(wildcards):
    return BUILDS[wildcards.build]

def get_fq(wildcards):
        return sequences.loc[(wildcards.unique_id), ["fq1", "fq2"]].dropna()

def get_rg(wildcards):
        return sequences.loc[(wildcards.unique_id), ["rg"]].dropna()["rg"]

def get_sample(wildcards):
        return sequences.loc[(wildcards.unique_id), ["Sample"]].dropna()["Sample"]

def getSampleGroupedBams(wildcards):
    return ["".join(wildcards.build + "/aligned_bam/"+wildcards.group+"/"+wildcards.group+"."+s+".clean.bam") for s in samplesInGroups[wildcards.group]]

def getSampleGroupedBais(wildcards):
    return ["".join(wildcards.build+"/aligned_bam/"+wildcards.group+"/"+wildcards.group+"."+s+".clean.bai") for s in samplesInGroups[wildcards.group]]

def get_ref(wildcards):
    return config[wildcards.build + "_ref"]

def get_mills(wildcards):
    return config[wildcards.build + "_mills"]

def get_indels(wildcards):
    return config[wildcards.build + "_indels"]

def get_dbsnp(wildcards):
    return config[wildcards.build + "_dbsnp"]

def get_omni(wildcards):
    return config[wildcards.build + "_omni"]

def get_kgp(wildcards):
    return config[wildcards.build + "_kgp"]

def get_hapmap(wildcards):
    return config[wildcards.build + "_hapmap"]
rule all:
    input:
        expand("qc/fastqc/{sample}.{unique_id}/", zip, sample = sequences["Sample"], unique_id = sequences["unique_id"]),
#        expand("qc/fastqc/{sample}.{unique_id}.R2_fastqc.zip", zip, sample = sequences["Sample"], unique_id = sequences["unique_id"]),
#        expand("b37/{chr}.exome.list", chr = CHRS_b37),
        expand(expand("{{build}}/stats/bwa/{sample}/{sample}.{unique_id}_{{build}}.idxstats.tsv", zip, sample = sequences['Sample'], unique_id =  sequences["unique_id"]), build = ["b37","b38"] ),
        expand("{build}/vqsr_vcf/{chr}_all_samples_genotyped_{build}.vqsr.snps.indels.vcf.gz", chr = CHRS, build = "b37"),
        expand(expand("{{build}}/aligned_bam/marked_dup/{sample}_{{build}}.clean.markdup.bam", zip, sample = sequences['Sample'], unique_id = sequences['unique_id']), build = "b38"),
        expand("b38/aligned_bam/bqsr/{sample}_b38.bqsr.bam", sample = sequences.Sample.unique()),
        expand("{build}/vqsr_vcf/{chr}_all_samples_genotyped_{build}.vqsr.snps.indels.vcf.gz", chr = CHRS, build = "b38"),



include: "rules/fastqc.smk"
include: "rules/bwa_align.smk"
include: "rules/sort_bam.smk"
include: "rules/stats_clean_merge.smk"
include: "rules/mark_dup.smk"
include: "rules/bqsr.smk"
include: "vcf.snake"
#include: "b38_vcf.snake"
