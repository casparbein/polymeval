rule run_rasusa:
    input:
        "raw_reads/{sample}.fastq.gz",
    output:
        "raw_reads/{sample}.downsampled.fastq.gz"
    threads:
        10
    resources:
        mem_mb=50000
    log:
        "logs/run_ranusa/{sample}.log"
    params:
        downsample_nts = downsample_nucs
    conda:
        "../envs/ranusa.yaml"
    shell:
        """
        rasusa \
        reads \
        -s 100 \
        -b {params.downsample_nts} \
        {input} \
        -o {output}
        """

rule run_rasusa_combinations:
    input:
        "raw_reads/{sample}.fastq.gz",
    output:
        "raw_reads/{sample}-{amount}.fastq.gz"
    threads:
        10
    resources:
        mem_mb=50000
    log:
        "logs/run_ranusa/{sample}-{amount}.log"
    params:
        downsample_nts = lambda wildcards: downsample_dict[wildcards.amount]
    conda:
        "../envs/ranusa.yaml"
    shell:
        """
        rasusa \
        reads \
        -s 100 \
        -b {params.downsample_nts} \
        {input} \
        -o {output}
        """