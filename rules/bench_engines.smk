happy_wrapper = f"{wrapper_versions['happy']}/bio/hap.py/hap.py"
tabix_wrapper_generic = f"{wrapper_versions['tabix']}/bio/tabix/index"

def q_vcf(wc):   return CALLERS[wc.caller]["vcf"].format(sample=wc.sample)
def q_tbi(wc):   return q_vcf(wc) + ".tbi"
def t_vcf(wc):   return BENCHMARKS[wc.truth]["vcf"]
def t_bed(wc):   return BENCHMARKS[wc.truth]["bed"]
def t_build(wc): return BENCHMARKS[wc.truth]["build"]


rule tabix_vcf:
    input:  
        "variants/{vcf}.vcf.gz"
    output: 
        "variants/{vcf}.vcf.gz.tbi"
    log:    
        "logs/tabix_vcf/{vcf}.log"
    params: "-p vcf"
    wrapper: 
        tabix_wrapper_generic

rule bench_happy:
    input:
        query         = q_vcf,
        truth         = t_vcf,
        truth_regions = t_bed,
        genome        = lambda wc: REF_PLAIN[t_build(wc)],
        genome_index  = lambda wc: REF_FAI[t_build(wc)],
    output:
        multiext("benchmarks/happy/{caller}/{truth}/{sample}/{sample}_results",
                 ".runinfo.json", ".vcf.gz", ".summary.csv", ".extended.csv",
                 ".metrics.json.gz", ".roc.all.csv.gz",
                 ".roc.Locations.INDEL.csv.gz", ".roc.Locations.INDEL.PASS.csv.gz",
                 ".roc.Locations.SNP.csv.gz", ".roc.tsv")
    params:
        engine = "vcfeval",
        prefix = lambda wc: f"benchmarks/happy/{wc.caller}/{wc.truth}/{wc.sample}/{wc.sample}_results",
        extra  = lambda wc: BENCHMARKS[wc.truth].get("happy_extra", "--verbose --pass-only"),
    log: 
        "logs/bench_happy/{caller}.{truth}.{sample}.log"
    threads: 4
    resources: 
        mem_mb = 200000
    wrapper: happy_wrapper


rule bench_truvari:
    input:
        comp  = q_vcf,
        index = q_tbi,
        base  = t_vcf,
        bed   = t_bed,
    output: 
        "benchmarks/truvari/{caller}/{truth}/{sample}/summary.json"
    params:
        out = "benchmarks/truvari/{caller}/{truth}/{sample}/",
        ref = lambda wc: REFERENCE[t_build(wc)],
        extra = lambda wc: BENCHMARKS[wc.truth].get("truvari_extra", "--sizemin 50"),
    log: 
        "logs/bench_truvari/{caller}.{truth}.{sample}.log"
    threads: 1
    resources: 
        mem_mb = 20000
    conda: 
        "../envs/truvari.yaml"
    shell:
        """
        rm -rf {params.out}
        truvari bench \
          -b {input.base} \
          -c {input.comp} \
          -o {params.out} \
          --includebed {input.bed} \
          --reference {params.ref} \
          {params.extra} 2> {log}
        """


rule bench_truvari_refine:
    input:
        summary = "benchmarks/truvari/{caller}/{truth}/{sample}/summary.json",
        comp    = q_vcf,
        index   = q_tbi,
    output: 
        "benchmarks/truvari/{caller}/{truth}/{sample}/refine.variant_summary.json"
    params:
        out = "benchmarks/truvari/{caller}/{truth}/{sample}/",
        ref = lambda wc: REFERENCE[t_build(wc)],
    log: 
        "logs/bench_truvari_refine/{caller}.{truth}.{sample}.log"
    threads: 1
    resources: 
        mem_mb = 200000
    conda: 
        "../envs/truvari.yaml"
    shell:
        """
        truvari refine \
          --use-original-vcfs \
          --reference {params.ref} \
          --buffer 0 \
          --coords O \
          --write-phab \
          {params.out} 2> {log}
        """


rule bench_aardvark:
    input:
        query       = q_vcf,
        query_index = q_tbi,
        truth       = t_vcf,
        regions     = t_bed,
    output: 
        "benchmarks/aardvark/{caller}/{truth}/{sample}/summary.tsv"
    params:
        out   = "benchmarks/aardvark/{caller}/{truth}/{sample}/",
        ref   = lambda wc: REFERENCE[t_build(wc)],
        gap   = lambda wc: 1000 if BENCHMARKS[wc.truth]["cls"] in ("sv", "tr") else 100,
        extra = lambda wc: "--enable-record-basepair-metrics"
                           if BENCHMARKS[wc.truth]["cls"] == "tr" else "",
    log: 
        "logs/bench_aardvark/{caller}.{truth}.{sample}.log"
    threads: 16
    resources: 
        mem_mb = 50000
    conda: 
        "../envs/aardvark.yaml"
    shell:
        """
        aardvark compare \
          --threads {threads} \
          --reference {params.ref} \
          --truth-vcf {input.truth} \
          --query-vcf {input.query} \
          --regions {input.regions} \
          -o {params.out} \
          --min-variant-gap {params.gap} \
          --compare-label {wildcards.caller} \
          {params.extra} 2> {log}
        """