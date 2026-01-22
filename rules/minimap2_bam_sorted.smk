## align raw reads with minimap2, sort, convert to bam and index
rule minimap2_bam_sorted:
    input:
        target=reference_seq,  # can be either genome index or genome fasta
        query="raw_reads/{sample}.fastq.gz",
    output:
        "alignments/{sample}.sorted.bam",
        idx = "alignments/{sample}.sorted.bam.csi"
    log:
        "logs/minimap2_bam_sorted/{sample}.log",
    params:
        extra=" -x map-hifi --MD --eqx -L ",  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard --threads maybe does not work, actually anything kills the downstream steps of the script
    resources:
        mem_mb = 200000
    threads: 50
    wrapper:
       # "v7.6.0/bio/minimap2/aligner"
        "https://github.com/snakemake/snakemake-wrappers/raw/release-please--branches--master/bio/minimap2/aligner"