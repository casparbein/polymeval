## Get read QC histogram with bbmap
rule nt_qc_hist:
    input:
        "raw_reads/{sample}.fastq.gz",
    output:
        out = temp("raw_reads/{sample}.tmp.fastq.gz"),
        qchist = "out/stats/{sample}.qchist.txt",
    log:
        "logs/nt_qc_hist/{sample}.log",
    params:
        command="reformat.sh",
    threads: 10
    resources:
        mem_mb=50000,
    wrapper:
        "v7.7.0/bio/bbtools"