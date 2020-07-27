rule sort_bam:
    input: 
        BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".bam"
    output: 
        bam=temp(BUILD + "/aligned_bam_"+ BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam"),
        bai=temp(BUILD + "/aligned_bam_" + BUILD + "/{sample}/{sample}.{unique_id}_" + BUILD + ".sorted.bam.bai")
    threads: 8
    shell:
        """
            samtools sort --threads {threads} -o {output.bam} {input} &&
            samtools index {output.bam}
        """

