## Find MDA artifacts and putative chimeras with breakinator
rule breakinator:
    input:
        bam = "alignments/{sample}.sorted.bam",
        ref = reference_seq
    output:
        summary_table = "out/files/{sample}_breakinator_summary.txt",
        detailed_breakinator = "out/files/{sample}_breakinator_detailed.txt",
    threads:
        3
    resources:
        mem_mb = 100000
    conda:
        "../envs/breakinator.yaml"
    log:
        "logs/breakinator/{sample}.log",
    shell:
        """
        breakinator \
        -i {input.bam} \
        -q 0 \
        -a 50 \
        --no-sym \
        -g {input.ref} \
        -c 50000 \
        -f 10000 \
        --tabular \
        -o {output.detailed_breakinator} > {output.summary_table}
        """