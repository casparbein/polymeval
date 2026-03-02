## Mapping against hg38 for SNVs/InDels
rule minimap2_bam_sorted:
    input:
        target=reference_seq,  #should be human hg38
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
    threads: 
        50
    wrapper:
       "v9.0.1/bio/minimap2/aligner"

## Mapping to hg37 for structural variants
rule minimap2_bam_sorted_hg37:
    input:
        target=reference_seq2, ## should by hg37 
        query="raw_reads/{sample}.fastq.gz",
    output:
        "alignments/{sample}.hs37d5.sorted.bam",
        idx = "alignments/{sample}.hs37d5.sorted.bam.csi"
    log:
        "logs/minimap2_bam_sorted_hg37/{sample}.log",
    params:
        extra=" -x map-hifi --MD --eqx -L ",  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard --threads maybe does not work, actually anything kills the downstream steps of the script
    resources:
        mem_mb = 200000
    threads: 
        50
    wrapper:
       "v9.0.1/bio/minimap2/aligner"