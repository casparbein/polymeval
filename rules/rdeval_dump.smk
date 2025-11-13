## Optional: dump read stats from deduplicated read sets with rdeval
rule rdeval_dump:
    input:
        "raw_reads/{sample}.fastq.gz",
    output:
        "out/stats/{sample}.rdeval_dump.tsv",
    threads:
        5
    resources:
        mem_mb = 100000
    log:
        "logs/rdeval_dump/{sample}.log",
    conda:
        "../envs/rdeval.yaml",
    shell:
        """
        rdeval \
        --sequence-report \
        -j {threads} \
        {input} > {output} \
        2> {log}
        """