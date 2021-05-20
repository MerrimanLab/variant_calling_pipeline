rule sample_coverage_stats:
    input:
        bam="{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.bam",
        bai="{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.bai"

    output:
        "{build}/stats/coverage/bqsr/{sample}_{build}.bqsr.coverage.tsv"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} samtools coverage {input.bam} > {output}"


rule fastq_coverage_stats:
    input:
        bam="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam",
        bai="{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam.bai"
    output:
        "{build}/stats/coverage/bwa/{sample}/{sample}.{unique_id}_{build}.coverage.tsv"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} samtools coverage {input.bam} > {output}"
