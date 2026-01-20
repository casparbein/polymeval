## Get average coverage/GC in windows of size 200 bp
rule pandepth_genome_small_window:
    input:
        ali = "alignments/{sample}.sorted.bam"
    output:
        "out/files/{sample}_pandepth_all.win.stat.gz"
    threads:
        16
    resources:
        mem_mb = 50000
    params:
        ref = reference_seq,
        out_prefix = 'out/files/{sample}_pandepth_all',
        pandepth_path = config["pandepth_path"]
    log:
        "logs/pandepth_genome_window/{sample}.log"
    shell:
        """
        {params.pandepth_path}/pandepth \
        -i {input.ali} \
        -d 0 \
        -t {threads} \
        -r {params.ref} \
        -c \
        -o {params.out_prefix} \
        -w 200 
        """