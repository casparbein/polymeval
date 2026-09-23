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