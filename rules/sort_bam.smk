rule sort_bam:
    input: 
         "{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.bam"
    output: 
        bam=temp("{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam"),
        bai=temp("{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.sorted.bam.bai")
    threads: 8
    shell:
        """
            samtools sort --threads {threads} -o {output.bam} {input} &&
            samtools index {output.bam}
        """

