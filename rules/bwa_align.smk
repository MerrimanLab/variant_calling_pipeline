rule align_fastq:
    input: 
        ref = get_ref,
        fq = get_fq     
    output: 
        temp("{build}/aligned_bam/{sample}/{sample}.{unique_id}_{build}.bam")
    log:
        "{build}/logs/bwamem/{sample}/{unique_id}.log"
    threads: 8
    params:
        nice = nice_cmd,
        RG= get_rg
    message: "Aligning unique_id: {wildcards.unique_id} using read group {params.RG}.\nInput is {input}.\nOutput file is {output}"
    shell:
        "{params.nice} bwa mem -t {threads} -R \"{params.RG}\" -M {input.ref} {input.fq} 2> {log} | samtools view -@ {threads} -Shb > {output}" 

