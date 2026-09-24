sniffles_wrapper = f"{wrapper_versions['sniffles']}/bio/sniffles"

## Call SVs with Sniffles2
rule run_sniffles2:
    input:
        samples = "alignments/{sample}.hs37d5.sorted.bam",
        ref = reference_seq2,
    output:
        vcf="variants/{sample}.sniffles.vcf.gz",
    log:
        "logs/run_sniffles2/{sample}.log"
    threads:
        4
    resources:
        mem_mb = 100000
    params:
        extra="", 
    wrapper:
        sniffles_wrapper

rule run_sniffles2_grch38:
    input:
        samples = "alignments/{sample}.sorted.bam",
        ref = reference_seq_gz,
    output: 
        vcf="variants/{sample}.sniffles.grch38.vcf.gz"
    log:   
        "logs/run_sniffles2_grch38/{sample}.log"
    threads: 4
    resources: 
        mem_mb = 100000
    params:
        extra="",  
    wrapper:
        sniffles_wrapper