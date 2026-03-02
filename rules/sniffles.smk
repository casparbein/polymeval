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
        --threads {threads}
        """