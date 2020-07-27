import pandas as pd

configfile: "config.yaml"

nice_cmd="nice -n 10"

sequences = pd.read_table(config["sequences_metadata"], dtype=str).set_index(["unique_id"], drop=False)
sequences["run_info"] = sequences.unique_id.str.replace(r'^([A-Za-z\d-]+_)', "")

samplesInGroups = {}
for sample in sequences["Sample"].unique():
  samplesInGroups[sample] = sequences.loc[(sequences['Sample'] == sample), ["unique_id"]]["unique_id"].to_list()

CHRS = list(range(1,23))

rule all:
    input:
        expand("b37/{chr}.exome.list", chr = CHRS),
        expand("b37/stats/bwa/{sample}/{sample}.{unique_id}_b37.idxstats.tsv", zip, sample = sequences['Sample'], unique_id =  sequences["unique_id"]),
        expand("b37/vqsr_vcf/{chr}_all_samples_genotyped_b37.vqsr.snps.indels.vcf.gz", chr = CHRS)





include: "b37.snake"


