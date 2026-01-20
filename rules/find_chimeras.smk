## Extract chimeras from reads
rule find_chimeras:
    input:
        bam = "alignments/{sample}.sorted.bam",
        fastq = "raw_reads/{sample}.fastq.gz"
    output:
        map_type = "out/files/{sample}_mappingType.txt",
        prim_hist = "out/files/{sample}_prim.hist.txt",
        supp_hist = "out/files/{sample}_suppl.hist.txt",
        rl_hist = "out/files/{sample}_rl.hist.txt",
    threads:
        1
    resources:
        mem_mb = 100000
    params:
        out_path = "out/files/{sample}"
    conda:
        "../envs/find_chimeras.yaml"
    log:
        "logs/find_chimeras/{sample}.log",
    script:
        "../scripts/find_chimeras.py"