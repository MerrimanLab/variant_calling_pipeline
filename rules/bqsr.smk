
rule base_recal:
    input:
        bam="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bam",
        bai="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bai",
        dbsnp = get_dbsnp,
        mills = get_mills,
        indels = get_indels,
        ref = get_ref
    output:
        bqsr = "{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.table"
    log:
        "{build}/logs/bqsr/{sample}.bqsr.log"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.dbsnp} --known-sites {input.mills} --known-sites {input.indels} -O {output} 2> {log}"


rule apply_bqsr:
    input:
        bam="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bam",
        bai="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bai",
        ref = get_ref,
        bqsr = "{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.table"
    output:
        bam="{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.bam",
        bai="{build}/aligned_bam/bqsr/{sample}_{build}.bqsr.bai"
    log:
        "{build}/logs/bqsr/{sample}.apply_bqsr.log"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.bqsr} -O {output.bam} 2> {log}"


