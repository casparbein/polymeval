meryl_wrapper =  f"{wrapper_versions['meryl']}/bio/meryl/count"

## Run meryl kmer counting for merqury 
rule meryl_count:
    input:
        fasta="raw_reads/{sample}.fastq.gz",
    output:
        temp(directory("meryl/{sample}_reads.meryl")),
    log:
        "logs/meryl_count/{sample}.log",
    params:
        command="count",
        extra="k=21",
    threads: 20
    resources:
        mem_mb=40000,
    wrapper:
        meryl_wrapper

## Run Merqury Assembly quality evaluation
rule run_merqury:
    input:
        db  = lambda wc: f"meryl/{asm_sample(wc.asm_id)}_reads.meryl",
        asm = "assemblies/{asm_id}.fa",
    output:
        "merqury/{asm_id}_slf/{asm_id}_slf.qv",
        temp(directory("merqury/{asm_id}_slf/{asm_id}.meryl")),
    threads:
        20
    resources:
        mem_mb=60000
    params:
        output_dir="merqury/{asm_id}_slf",
        out_prefix="{asm_id}_slf"
    log:
        "logs/run_merqury/{asm_id}.log",
    conda:
        "../envs/merqury.yaml"
    shell:
        """
        mkdir -p {params.output_dir};
        cd {params.output_dir};
        merqury.sh \
        ../../{input.db} \
        ../../{input.asm} \
        {params.out_prefix} \
        2> ../../{log};
        """