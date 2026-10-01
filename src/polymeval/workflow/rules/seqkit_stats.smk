seqkit_wrapper = f"{wrapper_versions['seqkit']}/bio/seqkit"

## General summary stats for (deduplicated) read sets
rule seqkit_stats:
    input:
        fastx=expand("raw_reads/{sample}.fastq.gz", sample = samples),
    output:
        stats="out/stats/seqkit_all.tsv",
    log:
        "logs/stats/seqkit.log",
    params:
        command="stats",
        extra="--all --tabular",
    threads: 10
    wrapper:
        seqkit_wrapper
