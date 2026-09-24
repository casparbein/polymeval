happy_wrapper = f"{wrapper_versions['happy']}/bio/hap.py/hap.py"
tabix_wrapper_generic = f"{wrapper_versions['tabix']}/bio/tabix/index"
bcftools_reheader_wrapper = f"{wrapper_versions['bcftools']}/bio/bcftools/reheader"
bcftools_view_wrapper     = f"{wrapper_versions['bcftools']}/bio/bcftools/view"

## For longcalld small vcf
def q_vcf(wc):
    c = CALLERS[wc.caller]
    key = "vcf_small" if BENCHMARKS[wc.truth]["cls"] == "small" and "vcf_small" in c else "vcf"
    return c[key].format(sample=wc.sample)

def q_tbi(wc):   return q_vcf(wc) + ".tbi"
def t_vcf(wc):   return f"benchmarks/truth/{wc.truth}.vcf.gz"
def t_tbi(wc):   return f"benchmarks/truth/{wc.truth}.vcf.gz.tbi"
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

## GIAB's SV Tier1 v0.6 defines BREAKSIMLENGTH twice. Reheader for aardvark to work
rule truth_header:
    input:  
        lambda wc: BENCHMARKS[wc.truth]["vcf"]
    output: 
        temp("benchmarks/truth/{truth}.hdr")
    conda:  
        "../envs/htslib.yaml"
    log:    
        "logs/truth_header/{truth}.log"
    shell:
        """
        bcftools view -h {input} | awk '
            /^##(INFO|FORMAT|FILTER)=<ID=/ {{
                tag = $0
                sub(/^##/, "", tag); sub(/=<ID=/, "\\t", tag); sub(/[,>].*$/, "", tag)
                if (tag in seen) next
                seen[tag] = 1
            }}
            {{ print }}' > {output} 2> {log}
        """

rule truth_vcf:
    input:
        lambda wc: BENCHMARKS[wc.truth]["vcf"],
        header = "benchmarks/truth/{truth}.hdr",
    output:
        "benchmarks/truth/{truth}.vcf.gz",
    params:
        extra = "",
    threads: 2
    log: 
        "logs/truth_vcf/{truth}.log"
    wrapper: 
        bcftools_reheader_wrapper

rule truth_tbi:
    input:   
        "benchmarks/truth/{truth}.vcf.gz"
    output:  
        "benchmarks/truth/{truth}.vcf.gz.tbi"
    params: 
        "-p vcf"
    log:     
        "logs/truth_tbi/{truth}.log"
    wrapper: 
        tabix_wrapper_generic

## longcalld emits small variants and SVs in one VCF, small-variant cells get a size-filtered copy.
rule small_only_vcf:
    input:   
        "variants/{sample}_longcalld.vcf.gz"
    output:  
        "variants/{sample}_longcalld.small.vcf.gz"
    params:  
        extra = "-e 'abs(ILEN)>=50'"
    log:     
        "logs/small_only_vcf/{sample}.log"
    wrapper: 
        bcftools_view_wrapper


## Happy benchmark
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

## Truvari
rule bench_truvari:
    input:
        comp  = q_vcf,
        index = q_tbi,
        base  = t_vcf,
        base_tbi = t_tbi,
        bed   = t_bed,
        ref      = lambda wc: REFERENCE[t_build(wc)],
        ref_fai  = lambda wc: REFERENCE[t_build(wc)] + ".fai",
        ref_gzi  = lambda wc: REFERENCE[t_build(wc)] + ".gzi",
    output: 
        "benchmarks/truvari/{caller}/{truth}/{sample}/summary.json"
    params:
        out = "benchmarks/truvari/{caller}/{truth}/{sample}/",
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
          --reference {input.ref} \
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
        truth_index = t_tbi,
        regions     = t_bed,
        ref      = lambda wc: REFERENCE[t_build(wc)],
        ref_fai  = lambda wc: REFERENCE[t_build(wc)] + ".fai",
        ref_gzi  = lambda wc: REFERENCE[t_build(wc)] + ".gzi",
    output: 
        "benchmarks/aardvark/{caller}/{truth}/{sample}/summary.tsv"
    params:
        out   = "benchmarks/aardvark/{caller}/{truth}/{sample}/",
        gap   = lambda wc: 1000 if BENCHMARKS[wc.truth]["cls"] in ("sv", "tr") else 100,
        extra = lambda wc: "--enable-record-basepair-metrics"
                           if BENCHMARKS[wc.truth]["cls"] == "tr" else "",
    log: 
        "logs/bench_aardvark/{caller}.{truth}.{sample}.log"
    threads: 4
    resources: 
        mem_mb = 50000
    conda: 
        "../envs/aardvark.yaml"
    shell:
        """
        aardvark compare \
          --threads {threads} \
          --reference {input.ref} \
          --truth-vcf {input.truth} \
          --query-vcf {input.query} \
          --regions {input.regions} \
          -o {params.out} \
          --min-variant-gap {params.gap} \
          --compare-label {wildcards.caller} \
          {params.extra} 2> {log}
        """

## Summarize all output tables
rule benchmark_table:
    input:
        [f for *_, f in BENCH_FILES]
    output:
        long   = "out/benchmarks/all_benchmarks.tsv",
        matrix = "out/benchmarks/benchmark_matrix.tsv",
    params:
        cells      = BENCH_FILES,
        benchmarks = BENCHMARKS,
    log:
        "logs/benchmark_table/collect.log"
    script:
        "../scripts/collect_benchmarks.py"