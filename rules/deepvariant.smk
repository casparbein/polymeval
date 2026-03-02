## Deepvariant small variant calling
rule run_deepvariant:
    input:
        reads = "alignments/{sample}.sorted.bam",
        reference = reference_seq,
    output:
        vcf="variants/{sample}.vcf.gz",
        gvcf = "variants/{sample}.gvcf.gz",
        #report="variants/{sample}.{reference_name}.visual_report.html"
    params:
        in_reads = "/input/{sample}.sorted.bam",
        in_ref = "/reference/{}".format(reference_seq.split('/')[-1]),
        out_vcf = "/output/{sample}.vcf.gz",
        out_gvcf = "/output/{sample}.gvcf.gz",
    container: "docker://google/deepvariant:1.10.0-beta",
    threads:
        32
    resources:
        mem_mb = 200000
    log:
        "logs/run_deepvariant/{sample}.log"
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref={params.in_ref} \
        --reads={params.in_reads} \
        --output_vcf={params.out_vcf} \
        --output_gvcf={params.out_gvcf} \
        --num_shards={threads} \
        --vcf_stats_report=true 
        """