def get_input_reads(wildcards):
    if config["gzipped"]:
        return "raw_reads/{sample}.dup.fastq.gz"
    else:
        return "raw_reads/{sample}.dup.fastq"

## Optional
rule mark_dups:
    input:
        reads=get_input_reads,#"raw_reads/{sample}.dup.fastq.gz",
    output:
        temp("raw_reads/{sample}.fastq"),
    threads:
        40
    resources:
        mem_mb = 100000
    params:
        dupreads="raw_reads/{sample}.indups.fastq", ## Change to fastqs
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
        40
    resources:
        mem_mb = 20000
    conda:
        "../envs/pigz.yaml",
    log:
        "logs/gzip_dedup/{sample}.log",
    shell:
        """
        pigz -p {threads} {input}
        """