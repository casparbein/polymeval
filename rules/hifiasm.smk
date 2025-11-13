#def hifieval_output(wilcards):
#    if 

## Polymerase-specific assembly with hifiasm
rule hifiasm:
    input:
        fasta="raw_reads/{sample}.fastq.gz",
    output:
        "assemblies/{sample}.p_ctg.gfa",
        "assemblies/{sample}.ec.fa" if config["hifieval"] else [], ## Error-corrected reads for hifieval
        temp("assemblies/{sample}.r_utg.gfa"),
        temp("assemblies/{sample}.p_utg.gfa"),
        temp("assemblies/{sample}.a_ctg.gfa"),
        temp("assemblies/{sample}.a_ctg.lowQ.bed"),
        temp("assemblies/{sample}.a_ctg.noseq.gfa"),
        temp("assemblies/{sample}.p_ctg.lowQ.bed"),
        temp("assemblies/{sample}.p_ctg.noseq.gfa"),
        temp("assemblies/{sample}.p_utg.lowQ.bed"),
        temp("assemblies/{sample}.p_utg.noseq.gfa"),
        temp("assemblies/{sample}.r_utg.lowQ.bed"),
        temp("assemblies/{sample}.r_utg.noseq.gfa"),
        temp("assemblies/{sample}.ec.bin"),
        temp("assemblies/{sample}.ovlp.reverse.bin"),
        temp("assemblies/{sample}.ovlp.source.bin"),
    log:
        "logs/hifiasm/{sample}.log",
    params:
        extra="--primary -l 3 --write-ec " if config["hifieval"] else " --primary -l 3  ",
    threads: 50
    resources:
        mem_mb=100000,
    wrapper:
        "v5.10.0/bio/hifiasm"


## Get Fasta file from gfa
rule get_fasta:
    input:
        "assemblies/{sample}.p_ctg.gfa"
    output:
        "assemblies/{sample}.fa"
    log:
        "logs/get_fasta/{sample}.log"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} > {output} 2> {log}
        """

## Get List of chromosome lengths for N(x) plotting
rule samtools_faidx:
    input:
        "assemblies/{sample}.fa",
    output:
        "assemblies/{sample}.fa.fai",
    log:
        "logs/faidx/{sample}.log",
    params:
        extra="",
    wrapper:
        "v7.6.0/bio/samtools/faidx"