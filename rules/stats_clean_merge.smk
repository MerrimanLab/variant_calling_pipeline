
rule idx_stats:
    input:
        bam=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam.bai"
    output:
        BUILD + "/stats/bwa/{sample}/{sample}.{unique_id}_" + BUILD + ".idxstats.tsv"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} samtools idxstats {input.bam} > {output}"


rule clean:
    input:
        bam=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam.bai"
    output:
        bam=temp(BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}.clean.bam"),
        bai=temp(BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}.clean.bai")
    log:
        BUILD + "/logs/cleansam/{sample}/{unique_id}_clean.log"
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
        bam=BUILD + "/aligned_bam_" + BUILD + "/{group}/{group}_" + BUILD + ".merged.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/{group}/{group}_" + BUILD + ".merged.bai"
    log:
        BUILD + "/logs/merge_bam/{group}/{group}_merge.log"
    params:
        nice = nice_cmd
    threads: 2
    shell:
        "{params.nice} gatk MergeSamFiles  $(echo {input.bam} | sed 's/^\| / -I /g')  -O {output.bam} --CREATE_INDEX true --ASSUME_SORTED true --USE_THREADING true 2> {log}"


