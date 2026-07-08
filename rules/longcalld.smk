tabix_wrapper = f"{wrapper_versions['tabix']}/bio/tabix/index"

## Run Tandem Repeat Genotyper with TRGT
rule run_longcalld:      
    input:
        ali = "alignments/{sample}.sorted.bam",
    output:
        "variants/{sample}_longcalld.vcf.gz",
    threads:
        16
    resources:
        mem_mb = 100000
    params:
        ref = reference_seq_gz,
    log:
        "logs/run_longcalld/{sample}.log"
    conda:
        "../envs/longcalld.yaml"
    shell:
        """
        longcallD \
        call \
        --hifi \
        -t {threads} \
        {params.ref} \
        {input.ali} \
        -o {output} \
        -Oz \
        2> {log}
        """

## Index VCF files
rule tabix_longcall:
    input:
        "variants/{sample}_longcalld.vcf.gz",
    output:
        "variants/{sample}_longcalld.vcf.gz.tbi",
    log:
        "logs/tabix_longcall/{sample}.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        tabix_wrapper