rule mark_dup:
    input:
        bam=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}_" + BUILD + ".merged.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}_" + BUILD + ".merged.bai"
    output:
        bam=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bam",
        bai=BUILD + "/aligned_bam_" + BUILD + "/marked_dup/{sample}_" + BUILD + ".clean.markdup.bai",
        metrics=BUILD + "/stats/mark_dups/{sample}.dup_metrics.txt"
    log:
        BUILD + "/logs/markdup/{sample}_markdup.log"
    params:
        nice = nice_cmd
    shell:
        """
            {params.nice} gatk MarkDuplicates -I {input.bam} -O {output.bam} -METRICS_FILE {output.metrics} -CREATE_INDEX true  -ASSUME_SORTED true 2> {log} 
        """

