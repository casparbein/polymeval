## Call SVs with Sniffles2
rule run_sniffles2:
    input:
        reads = "alignments/{sample}.hs37d5.sorted.bam",
        reference = reference_seq2,
    output:
        vcf="variants/{sample}.sniffles.vcf.gz",
    log:
        "logs/run_sniffles2/{sample}.log"
    conda:
        "../envs/sniffles.yaml"
    threads:
        20
    resources:
        mem_mb = 100000
    shell:
        """
        sniffles \
        --input {input.reads} \
        --reference {input.reference} \
        --vcf {output.vcf} \
        --threads {threads} \
        2> {log}
        """

rule run_sniffles2_grch38:
    input:
        reads = "alignments/{sample}.sorted.bam",
        reference = reference_seq_gz,
    output: 
        "variants/{sample}.sniffles.grch38.vcf.gz"
    conda: 
        "../envs/sniffles.yaml"
    log:   
        "logs/run_sniffles2_grch38/{sample}.log"
    threads: 20
    resources: 
        mem_mb = 100000
    shell:
        """
        sniffles \
        --input {input.reads} \
        --reference {input.reference} \
        --vcf {output} \
        --threads {threads} 2> {log}
        """