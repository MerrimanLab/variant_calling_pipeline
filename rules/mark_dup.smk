rule mark_dup:
    input:
        bam="{build}/aligned_bam/{sample}/{sample}_{build}.merged.bam",
        bai="{build}/aligned_bam/{sample}/{sample}_{build}.merged.bai"
    output:
        bam="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bam",
        bai="{build}/aligned_bam/marked_dup/{sample}_{build}.clean.markdup.bai",
        metrics="{build}/stats/mark_dups/{sample}.dup_metrics.txt"
    log:
        "{build}/logs/markdup/{sample}_markdup.log"
    params:
        nice = nice_cmd,
        extra ="-CREATE_INDEX true  -ASSUME_SORTED true" 
    shell:
        """
            {params.nice} gatk MarkDuplicates -I {input.bam} -O {output.bam} -METRICS_FILE {output.metrics} {params.extra} 2> {log} 
        """

