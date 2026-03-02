## Run Truvari for SV benchmarking
rule run_truvari_sv:      
    input:
        comp = "variants/{sample}.sniffles.vcf.gz"
    output:
        "benchmarks/{sample}_truvari_sv/summary.json"
    threads:
        1
    resources:
        mem_mb = 100000
    params:
        base = sv_ground_vcf,
        bed = sv_ground_bed,
        out_prefix = 'benchmarks/{sample}_truvari_sv/',
    log:
        "logs/run_truvari_sv/{sample}.log"
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.out_prefix}; \
        truvari \
        bench \
        -b {params.base}\
        -c {input.comp} \
        -o {params.out_prefix} \
        --includebed {params.bed} \
        --passonly \
        -r 1000 \
        -p 0.0 \
        --sizemin 50
        """