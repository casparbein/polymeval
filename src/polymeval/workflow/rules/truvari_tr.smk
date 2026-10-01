## Benchmarking of tandem repeats with Truvari
## Setup based on: https://github.com/ACEnglish/truvari/issues/304
rule run_truvari_tr_bench:      
    input:
        comp = "variants/{sample}_trgt.sorted.norm.vcf.gz",
        index = "variants/{sample}_trgt.sorted.norm.vcf.gz.tbi",
    output:
        "benchmarks/{sample}_truvari_tr/summary.json"
    threads:
        1
    resources:
        mem_mb = 100000
    params:
        base = tandem_repeats_truvari, #HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.vcf
        base_bed = tandem_repeats_truvari_bed, #HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz
        out_prefix = 'benchmarks/{sample}_truvari_tr',
    log:
        "logs/run_truvari_tr_bench/{sample}.log"
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.out_prefix};
        truvari \
        bench \
        -b {params.base}\
        -c {input.comp} \
        -o {params.out_prefix} \
        --includebed {params.base_bed} \
        --sizemin 5 \
        --pick ac \
        2> {log}
        """

## Because of variant representation differences, phab is the right way to go
rule run_truvari_tr_refine:      
    input:
        comp = "variants/{sample}_trgt.sorted.norm.vcf.gz",
        index = "variants/{sample}_trgt.sorted.norm.vcf.gz.tbi",
        input_sum = "benchmarks/{sample}_truvari_tr/summary.json",
    output:
        "benchmarks/{sample}_truvari_tr/refine.variant_summary.json"
    threads:
        1
    resources:
        mem_mb = 100000
    params:
        reference = reference_seq_gz,
        bench_dir = "benchmarks/{sample}_truvari_tr/"
    log:
        "logs/run_truvari_tr/{sample}.log"
    conda:
        "../envs/truvari.yaml"
    shell:
        """
        truvari \
        refine \
        --use-original-vcfs \
        --reference {params.reference} \
        --buffer 0 \
        --coords O {params.bench_dir} \
        --write-phab \
        2> {log}
        """