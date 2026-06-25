rule prepare_bed_for_trgt:
    input:
        tandem_repeats_trgt, ##HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz
    output:
        "variants/tandem_repeats_for_trgt.bed",
    log:
        "logs/prepare_bed_for_trgt/out.log"
    shell:
        """
        gunzip -c {input} | cut -f1-3 | awk '{{i++; print $1 \"\t\" $2 \"\t\" $3 \"\tID=\" i \";MOTIFS=N;STRUC=(N)n\"}}' > {output}
        """


## Run Tandem Repeat Genotyper with TRGT
rule run_trgt:      
    input:
        ali = "alignments/{sample}.sorted.bam",
        tr_bed = "variants/tandem_repeats_for_trgt.bed", #tandem_repeats,
    output:
        "variants/{sample}_trgt.vcf.gz",
    threads:
        16
    resources:
        mem_mb = 100000
    params:
        ref = reference_seq,
        out_prefix = 'variants/{sample}_trgt',
    log:
        "logs/run_trgt/{sample}.log"
    conda:
        "../envs/trgt.yaml"
    shell:
        """
        trgt \
        genotype \
        --reads {input.ali} \
        --repeats {input.tr_bed} \
        --threads {threads} \
        --genome {params.ref} \
        --output-prefix {params.out_prefix} \
        2> {log}
        """

#bcftools norm -N -m-any repliQa_trgt.vcf.gz -O z -o comp.vcf.gz

rule bcftools_sort:
    input:
        "variants/{sample}_trgt.vcf.gz",
    output:
        "variants/{sample}_trgt.sorted.vcf.gz",
    log:
        "logs/bcftools/sort/{sample}.log",
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        #extras="-Oz",
    resources:
        mem_mb=8000,
    wrapper:
        f"{wrapper_versions['bcftools']}/bio/bcftools/sort"

rule norm_vcf:
    input:
        "variants/{sample}_trgt.sorted.vcf.gz",
        ref=reference_seq  # optional reference (will be translated into the -f option)
    output:
        "variants/{sample}_trgt.sorted.norm.vcf.gz",  # can also be .bcf, corresponding --output-type parameter is inferred automatically
    log:
        "logs/bcftools/norm/{sample}.norm.log",
    params:
        extra=" -N -m-any ",  # optional
        #uncompressed_bcf=False,
    wrapper:
        f"{wrapper_versions['bcftools']}/bio/bcftools/norm"

## Index VCF files
rule tabix:
    input:
        "variants/{sample}_trgt.sorted.norm.vcf.gz",
    output:
        "variants/{sample}_trgt.sorted.norm.vcf.gz.tbi",
    log:
        "logs/tabix/{sample}.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    wrapper:
        f"{wrapper_versions['tabix']}/bio/tabix/index"