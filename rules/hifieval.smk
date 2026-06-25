## Approximation of read error stats with hifieval
## Align raw reads
rule hifieval_align_raw:
    input:
        target="assemblies/{sample}.fa",  # can be either genome index or genome fasta
        query="raw_reads/{sample}.fastq.gz",
    output:
        temp("alignments/{sample}.raw.paf"),
    log:
        "logs/hifieval_align_raw/{sample}.log",
    params:
        extra="-cx map-hifi --secondary=no --paf-no-hit --cs", 
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 20
    resources:
        mem_mb = 100000
    wrapper:
       f"{wrapper_versions['minimap']}/bio/minimap2/aligner"

## Aligned error corrected reads
rule hifieval_align_ec:
    input:
        target="assemblies/{sample}.fa",  # can be either genome index or genome fasta
        query="assemblies/{sample}.ec.fa",
    output:
        temp("alignments/{sample}.ec.paf"),
    log:
        "logs/hifieval_align_ec/{sample}.log",
    params:
        extra="-cx map-hifi --secondary=no --paf-no-hit --cs", 
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 20
    resources:
        mem_mb = 100000
    wrapper:
        f"{wrapper_versions['minimap']}/bio/minimap2/aligner"

## Normally optional: If read file contains empty reads (Revio demo from PacBio does), hifieval will choke on this:
rule remove_empty:
    input:
        ec = "alignments/{sample}.ec.paf",
        raw = "alignments/{sample}.raw.paf",
    output:
        ec = temp("alignments/{sample}.ec.clean.paf"),
        raw = temp("alignments/{sample}.raw.clean.paf"),
    log: 
        "logs/remove_empty/{sample}.log",
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        awk '$5 != "*" {{print}}' {input.ec} > {output.ec};
        awk '$5 != "*" {{print}}' {input.raw} > {output.raw};
        """

## Run Hifieval
rule hifieval_compare:
    input:
        raw ="alignments/{sample}.raw.clean.paf", 
        ec ="alignments/{sample}.ec.clean.paf", 
    output:
        metric = "hifieval/{sample}.metric.eval.tsv",
        rdl_eval = "hifieval/{sample}.rdlvl.eval.tsv",
        summary = temp("hifieval/{sample}.summary.tsv"),
    params:
        out_base = "hifieval/{sample}",
    threads:
        1
    resources:
        mem_mb = 50000,
    log:
        "logs/hifieval_compare/{sample}.log",
    conda:
        "../envs/hifieval.yaml",
    shell:
        """
        hifieval.py \
        -o {params.out_base} \
        -r {input.raw} \
        -c {input.ec} \
        2> {log}
        """