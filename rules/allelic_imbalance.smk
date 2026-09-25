Alib_TRUTH_VCF = v5_smvar_vcf
Alib_TRUTH_BED = v5_smvar_bed

## Extract biallelic sites from truth
rule alib_truth_hets:
    input:
        vcf = Alib_TRUTH_VCF,
        tbi = Alib_TRUTH_VCF + ".tbi",
        bed = Alib_TRUTH_BED,
        ref = reference_seq,
    output:
        vcf = "allelic_imbalance/truth_het.vcf.gz",
        tbi = "allelic_imbalance/truth_het.vcf.gz.tbi",
    log: 
        "logs/alib_truth_hets/truth.log"
    conda: 
        "../envs/htslib.yaml"
    shell:
        """
        set -o pipefail
        bcftools view -R {input.bed} -m2 -M2 -g het {input.vcf} 2> {log} \
        | bcftools norm -f {input.ref} -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

## Write these into tsv, with allele depth for ref/alt
rule alib_truth_sites:
    input:
        vcf = "allelic_imbalance/truth_het.vcf.gz",
        tbi = "allelic_imbalance/truth_het.vcf.gz.tbi",
    output: 
        "allelic_imbalance/truth_het_sites.tsv"
    log: 
        "logs/alib_truth_sites/truth.log"
    conda: 
        "../envs/htslib.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} > {output} 2> {log}"

## Truth GC for het sites
rule alib_site_gc:
    input:
        sites = "allelic_imbalance/truth_het_sites.tsv",
        ref   = reference_seq, 
    output: 
        "allelic_imbalance/truth_het_sites.gc.tsv"
    params:
        flank = config.get("ai_gc_flank", 100),
    log: 
        "logs/alib_site_gc/truth.log"
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        set -o pipefail
        awk -v f={params.flank} 'BEGIN{{OFS="\t"}} \
              {{s = $2 - 1 - f; if (s < 0) s = 0; print $1, s, $2 + f, $1"_"$2}}' \
            {input.sites} | sort -k1,1 -k2,2n > {output}.bed 2> {log}
        bedtools nuc -fi {input.ref} -bed {output}.bed 2>> {log} \
          | tail -n +2 | cut -f4,6 > {output}
        rm -f {output}.bed
        """

## Get calls for query deepvariant vcf files at truth positions
rule alib_query_at_hets:
    input:
        query = "variants/{sample}_longcalld.small.vcf.gz",
        qidx  = "variants/{sample}_longcalld.small.vcf.gz.tbi",
        sites = "allelic_imbalance/truth_het.vcf.gz",
        sidx  = "allelic_imbalance/truth_het.vcf.gz.tbi",
        ref   = reference_seq,
    output: 
        "allelic_imbalance/{sample}.query_at_hets.tsv"
    log: 
        "logs/alib_query_at_hets/{sample}.log"
    conda: 
        "../envs/htslib.yaml"
    shell:
        """
        set -o pipefail
        printf 'chrom\tpos\tref\talt\tfilter\tgt\tdp\tad\tgq\n' > {output}
        bcftools norm -f {input.ref} -m -any {input.query} 2>> {log} \
        | bcftools view -T {input.sites} 2>> {log} \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%AD\t%GQ]\n' \
          >> {output} 2>> {log}
        """

## Coverage from samtools depth -> maybe just use het calls
rule alib_site_depth:
    input:
        bam   = "alignments/{sample}.sorted.bam",
        csi   = "alignments/{sample}.sorted.bam.csi",
        sites = "allelic_imbalance/truth_het_sites.tsv",
    output: 
        "allelic_imbalance/{sample}.site_depth.tsv"
    params:
        min_bq = config.get("ai_min_bq", 20),
        min_mq = config.get("ai_min_mq", 1),
    log: 
        "logs/alib_site_depth/{sample}.log"
    conda: 
        "../envs/htslib.yaml"
    shell:
        """
        set -o pipefail
        awk 'BEGIN{{OFS="\t"}} {{print $1, $2 - 1, $2}}' {input.sites} \
          | sort -k1,1 -k2,2n > {output}.bed
        samtools depth -a -b {output}.bed -q {params.min_bq} -Q {params.min_mq} \
            {input.bam} > {output} 2> {log}
        rm -f {output}.bed
        """

rule alib_analysis:
    input:
        query  = expand("allelic_imbalance/{sample}.query_at_hets.tsv", sample=samples),
        depth  = expand("allelic_imbalance/{sample}.site_depth.tsv",    sample=samples),
        sites  = "allelic_imbalance/truth_het_sites.tsv",
        gc     = "allelic_imbalance/truth_het_sites.gc.tsv",
    output:
        vaf     = "allelic_imbalance/alib_vaf_summary.tsv",
        dropout = "allelic_imbalance/alib_dropout.tsv",
        sig     = "allelic_imbalance/alib_significant_sites.tsv",
        plots   = "allelic_imbalance/alib_plots.pdf",
    params:
        sample_names = samples,
        colors       = config["colors"],
        min_dp       = config.get("ai_min_dp", 1),
        min_gq       = config.get("ai_min_gq", 1),
    resources: 
        mem_mb = 40000
    log: "logs/alib_analysis/analysis.log"
    conda: 
        "../envs/alib_stats.yaml"
    script: 
        "../scripts/allelic_imbalance.R"