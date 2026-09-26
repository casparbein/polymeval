## Run Tandem Repeat Genotyper with TRGT
rule run_longcalld:      
    input:
        ali = "alignments/{sample}.sorted.bam",
    output:
        "variants/{sample}_longcalld.vcf.gz",
    threads:
        16
    resources:
        mem_mb = 100000
    params:
        ref = reference_seq_gz,
    log:
        "logs/run_longcalld/{sample}.log"
    conda:
        "../envs/longcalld.yaml"
    shell:
        """
        longcallD \
        call \
        --hifi \
        -t {threads} \
        {params.ref} \
        {input.ali} \
        -o {output} \
        -Oz \
        2> {log}
        """

## longcalld emits small variants and SVs in one VCF, small-variant cells get a size-filtered copy.
rule small_only_vcf:
    input:   
        "variants/{sample}_longcalld.vcf.gz"
    output:  
        "variants/{sample}_longcalld.small.vcf.gz"
    params:  
        extra = "-i 'abs(ILEN)<50'"
    log:     
        "logs/small_only_vcf/{sample}.log"
    wrapper: 
        bcftools_view_wrapper

## Check pass filtering here later
rule small_only_vcf:
    input:   "variants/{sample}_longcalld.pass.vcf.gz"
    output:  "variants/{sample}_longcalld.pass.small.vcf.gz"
    params:  extra = "-i 'abs(ILEN)<50'"
    log:     "logs/small_only_vcf/{sample}.log"
    wrapper: bcftools_view_wrapper

rule sv_only_vcf:
    input:   
        "variants/{sample}_longcalld.pass.vcf.gz"
    output:  
        "variants/{sample}_longcalld.pass.sv.vcf.gz"
    params:  
        extra = "-e 'abs(ILEN)<50'"
    log:     
        "logs/sv_only_vcf/{sample}.log"
    wrapper: 
        bcftools_view_wrapper