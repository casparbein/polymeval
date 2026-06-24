## Get average coverage/GC in windows of size 200 bp
rule pandepth_genome_small_window:
    input:
        ali = "alignments/{sample}.sorted.bam",
        ref = reference_seq,
    output:
        "out/files/{sample}_pandepth_all.win.stat.gz"
    threads:
        16
    resources:
        mem_mb = 50000
    params:
        out_prefix = 'out/files/{sample}_pandepth_all',
        pandepth_path = config["pandepth_path"] 
    log:
        "logs/pandepth_genome_window/{sample}.log"
    #conda:
    #    "../envs/pandepth.yaml"
    container:
        "docker://ghcr.io/casparbein/pandepth:2.26"
    shell:
        """
        pandepth \
        -i {input.ali} \
        -d 0 \
        -t {threads} \
        -r {input.ref} \
        -c \
        -o {params.out_prefix} \
        -w 200 
        """