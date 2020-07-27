rule align_fastq:
    input: 
        ref = ref,
        fq = get_fq     
    output: 
        temp(BUILD + "/aligned_bam_" + BUILD +"/{sample}/{sample}.{unique_id}_" + BUILD + ".bam")
    log:
        BUILD + "/logs/bwamem/{sample}/{unique_id}.log"
    threads: 8
    params:
        nice = nice_cmd,
        RG= get_rg
    message: "Aligning unique_id: {wildcards.unique_id} using read group {params.RG}.\nInput is {input}.\nOutput file is {output}"
    shell:
        "{params.nice} bwa mem -t {threads} -p -R \"{params.RG}\" -M {input.ref} {input.fq} 2> {log} | samtools view -@ {threads} -Shb > {output}" 

