rule truvari_cmrg_sv:
    input:
        comp  = "variants/{sample}.sniffles.grch38.vcf.gz",
        index = "variants/{sample}.sniffles.grch38.vcf.gz.tbi",
    output: "benchmarks/{sample}_truvari_cmrg_sv/summary.json"
    params:
        base       = cmrg_sv_vcf,
        bed        = cmrg_sv_bed,
        ref        = reference_seq_gz,
        out_prefix = "benchmarks/{sample}_truvari_cmrg_sv/",
    conda: 
        "../envs/truvari.yaml"
    log: 
        "logs/truvari_cmrg_sv/{sample}.log"
    shell:
        """
        rm -rf {params.out_prefix};
        truvari \
        bench \
        -b {params.base} \
        -c {input.comp} \
        -o {params.out_prefix} \
        --includebed {params.bed} \
        --reference {params.ref} \
        --dup-to-ins \
        --sizemin 50 2> {log}
        """