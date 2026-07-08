## Run Tandem Repeat Genotyper with TRGT
rule run_aardvark:      
    input:
        truth = GIAB_joint_variants,
        truth_regions = GIAB_stvar_regions_bed,
        query = "variants/{sample}_longcalld.vcf.gz",
        query_index = "variants/{sample}_longcalld.vcf.gz.tbi",
    output:
        "benchmarks/aardvark/{sample}_aardvark/summary.tsv",
    threads:
        16
    resources:
        mem_mb = 100000
    params:
        ref = reference_seq_gz,
        out_prefix = 'benchmarks/aardvark/{sample}_aardvark/',
    log:
        "logs/run_aardvark/{sample}.log"
    conda:
        "../envs/aardvark.yaml"
    shell:
        """
        aardvark \
        compare  \
        --threads {threads} \
        --reference {params.ref} \
        --truth-vcf {input.truth} \
        --query-vcf {input.query} \
        --regions {input.truth_regions} \
        -o {params.out_prefix} \
        --min-variant-gap 1000 \
        --enable-record-basepair-metrics \
        2> {log}
        """