rule fastqc:
    input:
        fq = get_fq
    output:
        directory("qc/fastqc/{sample}.{unique_id}/")
    log:
        "logs/fastqc/{sample}.{unique_id}.log"
    benchmark:
        "benchmarks/fastqc/{sample}.{unique_id}.fastqc"
    threads: 8
    params:
        nice = nice_cmd
    message:
        "Undertaking quality control checks on raw sequence data"
    shell:
        "{params.nice} fastqc {input} --outdir {output} --threads {threads} 2> {log}"
