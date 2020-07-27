
rule base_recal:
    input:
        bam=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bai",
        dbsnp = dbsnp,
        mills = mills,
        indels = indels,
        ref = ref
    output:
        bqsr = BUILD + "/aligned_bam_" + BUILD + "/bqsr/{sample}_" + BUILD + ".bqsr.table"
    log:
        BUILD + "/logs/bqsr/{sample}.bqsr.log"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.dbsnp} --known-sites {input.mills} --known-sites {input.indels} -O {output} 2> {log}"


rule apply_bqsr:
    input:
        bam=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bai",
        ref = ref,
        bqsr = BUILD + "/aligned_bam_" + BUILD + "/bqsr/{sample}_" + BUILD + ".bqsr.table"
    output:
        bam=BUILD + "/aligned_bam_" + BUILD + "/bqsr/{sample}_" + BUILD + ".bqsr.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/bqsr/{sample}_" + BUILD + ".bqsr.bai"
    log:
        BUILD + "/logs/bqsr/{sample}.apply_bqsr.log"
    params:
        nice = nice_cmd
    shell:
        "{params.nice} gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.bqsr} -O {output.bam} 2> {log}"


