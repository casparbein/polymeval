## Run meryl kmer counting for merqury 
rule meryl_count:
    input:
        fasta="raw_reads/{sample}.fastq.gz",
    output:
        directory("meryl/{sample}_reads.meryl"),
    log:
        "logs/meryl_count/{sample}.log",
    params:
        command="count",
        extra="k=21",
    threads: 20
    resources:
        mem_mb=40000,
    wrapper:
        "v7.6.0/bio/meryl/count"

## Run Merqury Assembly quality evaluation
rule run_merqury:
    input:
        db="meryl/{sample}_reads.meryl",
        asm="assemblies/{sample}.fa",
    output:
        "merqury/{sample}_slf/{sample}_slf.qv",
    threads:
        20
    resources:
        mem_mb=60000
    params:
        output_dir="merqury/{sample}_slf",
        out_prefix="{sample}_slf"
    log:
        "logs/run_merqury/{sample}.log",
    conda:
        "../envs/merqury.yaml"
    shell:
        """
        mkdir -p {params.output_dir};
        cd {params.output_dir};
        merqury.sh \
        ../../{input.db} \
        ../../{input.asm} \
        {params.out_prefix};
        """