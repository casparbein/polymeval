## Optional
rule mark_dups:
    input:
        reads="raw_reads/{sample}.dup.fastq.gz",
    output:
        temp("raw_reads/{sample}.fastq"),
    threads:
        40
    resources:
        mem_mb = 100000
    params:
        dupreads="raw_reads/{sample}.indups.fastq.gz", ## Change to fastqs
    conda:
        "../envs/pbmarkdup.yaml",
    log:
        "logs/mark_dups/{sample}.log",
    shell:
        """
        pbmarkdup \
        --dup-file {params.dupreads} \
        -j {threads} \
        {input.reads} \
        --log-level INFO \
        {output} \
        2> {log} 
        """

## Gzip deduped reads
rule gzip_dedup:
    input:
        "raw_reads/{sample}.fastq",
    output:
        "raw_reads/{sample}.fastq.gz",
    threads:
        1
    resources:
        mem_mb = 20000
    log:
        "logs/gzip_dedup/{sample}.log",
    shell:
        """
        gzip {input}
        """