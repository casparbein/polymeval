happy_wrapper = f"{wrapper_versions['happy']}/bio/hap.py/hap.py"

## Running happy on deepvariant data
rule benchmark_variants:
    input:
        truth=vcf_benchmark,
        query="variants/{sample}.vcf.gz",
        truth_regions=vcf_truth_bed,
        #strats="stratifications.tsv",
        #strat_dir="strats_dir",
        genome=reference_seq,
        genome_index=reference_seq_idx,
    output:
        multiext("benchmarks/happy/{sample}_happy/{sample}_results",".runinfo.json",".vcf.gz",".summary.csv",
                ".extended.csv",".metrics.json.gz",".roc.all.csv.gz",
                ".roc.Locations.INDEL.csv.gz",".roc.Locations.INDEL.PASS.csv.gz",
                ".roc.Locations.SNP.csv.gz",".roc.tsv")
    params:
        engine="vcfeval",
        prefix=lambda wc: f"benchmarks/happy/{wc.sample}_happy/{wc.sample}_results",
        #prefix=lambda wc, input, output: output[0].split('.')[0],
        ## parameters such as -L to left-align variants
        extra="--verbose --pass-only"
    log: 
        "logs/benchmark_variants/{sample}.log"
    threads: 20
    resources:
        mem_mb = 600000
    wrapper: 
        happy_wrapper