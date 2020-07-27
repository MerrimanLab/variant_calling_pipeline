import pandas as pd

configfile: "config.yaml"

CHRS = list(range(1,23))



nice_cmd="nice -n 10"

sequences = pd.read_table(config["sequences_metadata"], dtype=str).set_index(["unique_id"], drop=False)
sequences["run_info"] = sequences.unique_id.str.replace(r'^([A-Za-z\d-]+_)', "")

samplesInGroups = {}
for sample in sequences["Sample"].unique():
  samplesInGroups[sample] = sequences.loc[(sequences['Sample'] == sample), ["unique_id"]]["unique_id"].to_list()


def get_fq(wildcards):
        return sequences.loc[(wildcards.unique_id), ["fq1", "fq2"]].dropna()

def get_rg(wildcards):
        return sequences.loc[(wildcards.unique_id), ["rg"]].dropna()["rg"]

def get_sample(wildcards):
        return sequences.loc[(wildcards.unique_id), ["Sample"]].dropna()["Sample"]


def getSampleGroupedBams(wildcards):
    return ["".join("b37/aligned_bam_b37/"+wildcards.group+"/"+wildcards.group+"."+s+".clean.bam") for s in samplesInGroups[wildcards.group]]
    
def getSampleGroupedBais(wildcards):
    return ["".join("b37/aligned_bam_b37/"+wildcards.group+"/"+wildcards.group+"."+s+".clean.bai") for s in samplesInGroups[wildcards.group]]




rule all:
    input:
        expand("b37/{chr}.exome.list", chr = CHRS),
        expand("b37/stats/bwa/{sample}/{sample}.{unique_id}_b37.idxstats.tsv", zip, sample = sequences['Sample'], unique_id =  sequences["unique_id"]),
        expand("b37/vqsr_vcf/{chr}_all_samples_genotyped_b37.vqsr.snps.indels.vcf.gz", chr = CHRS)#,
#        expand("b38/aligned_bam_b38/bqsr/{sample}_b38.bqsr.bam", sample = sequences.Sample.unique())





include: "b37_align.snake"
include: "b37_vcf.snake"
#include: "b38_align.snake"
#include: "b38_vcf.snake"


