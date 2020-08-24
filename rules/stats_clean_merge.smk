
rule idx_stats:
    input:
        bam="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam",
        bai="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam.bai"
    output:
        "{build}/stats/bwa/{sample}/{sample}.{unique_id}_{build}.idxstats.tsv"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} samtools idxstats {input.bam} > {output}"


rule clean:
    input:
        bam="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam",
        bai="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam.bai"
    output:
        bam=temp("{build}/aligned_bam/{sample}/{sample}.{unique_id}.clean.bam"),
        bai=temp("{build}/aligned_bam/{sample}/{sample}.{unique_id}.clean.bai")
    log:
        "{build}/logs/cleansam/{sample}/{unique_id}_clean.log"
    wildcard_constraints:
        sample="[A-Za-z\d-]+"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} gatk CleanSam -I {input.bam} -O {output.bam} -CREATE_INDEX true 2> {log}"

rule merge_bams:
    input:
        bam= getSampleGroupedBams,
        bai= getSampleGroupedBais
    output:
        bam= "{build}/aligned_bam/{group}/{group}_{build}.merged.bam",
        bai="{build}/aligned_bam/{group}/{group}_{build}.merged.bai"
    log:
        "{build}/logs/merged_bam/{group}/{group}_merge.log"
    params:
        nice = nice_cmd,
        extra ="--CREATE_INDEX true --ASSUME_SORTED true --USE_THREADING true" 
    threads: 2
    shell:
        "{params.nice} gatk MergeSamFiles  $(echo {input.bam} | sed 's/^\| / -I /g')  -O {output.bam} {params.extra}  2> {log}"


