## Rule to get KMC kmer dumps
## -hc was removed for now
rule kmc_count:
    input:
        "raw_reads/{sample}.fastq.gz"
    output:
        suf = "kmc/{sample}.res.kmc_suf",
        pre = "kmc/{sample}.res.kmc_pre"
    threads:
        10
    resources:
        mem_mb  = 40000
    params:
        tmp_dir = "kmc/{sample}_tmp_dir/",
        out_base = "kmc/{sample}.res",
        memory=40
    log:
        "logs/kmc_count/{sample}.log",
    envmodules:
        "kmc/3.2.4"
    shell:
        """
        mkdir -p {params.tmp_dir};
        kmc \
        -v \
        -k25 \
        -ci1 \
        -cs5000 \
        -t{threads} \
        -m{params.memory} \
        {input} \
        {params.out_base} \
        {params.tmp_dir} \
        --opt-out-size \
        2> {log}
        """

## Dump histogram into hist file (for example for genomescope)
rule kmc_dump_hist:
    input:
        suf = "kmc/{sample}.res.kmc_suf",
        pre = "kmc/{sample}.res.kmc_pre"
    output:
        "out/stats/{sample}.kmc_hist.txt"
    threads:
        10
    resources:
        mem_mb  = 40000
    params:
        db = "kmc/{sample}.res",
        memory=40
    log:
        "logs/kmc_dump_hist/{sample}.log",
    envmodules:
        "kmc/3.2.4"
    shell:
        """
        kmc_tools \
        transform \
        {params.db} \
        -ci0 \
        -cx5000 \
        histogram \
        {output} \
        2> {log}
        """
## Build unions/intersections/differences between K-mer datasets
rule kmc_union:
    input:
        kmc = expand("kmc/{sample}.res.kmc_suf", sample = samples),
    output:
        script = "kmc/union_script",
        kmc_out_pre = "kmc/kmc_union.kmc_pre",
        kmc_out_suf = "kmc/kmc_union.kmc_suf"
    threads:
        10
    resources:
        mem_mb  = 40000
    params:
        memory=40
    log:
        "logs/kmc_union/log",
    envmodules:
        "kmc/3.2.4"
    shell:
        """
        python3 /beegfs/home/bbein/polymerase_snakemake_wrapper/prepare_kmc.py \
        -i {input.kmc} \
        -o {output.script} \
        -u;
        kmc_tools \
        complex \
        {output.script} \
        2> {log}
        """

## Build intersection
rule kmc_isec:
    input:
        kmc = expand("kmc/{sample}.res.kmc_suf", sample = samples),
    output:
        script = "kmc/isec_script",
        kmc_out_pre = "kmc/kmc_intersection.kmc_pre",
        kmc_out_suf = "kmc/kmc_intersection.kmc_suf"
    threads:
        10
    resources:
        mem_mb  = 40000
    params:
        memory=40
    log:
        "logs/kmc_isec/log",
    envmodules:
        "kmc/3.2.4"
    shell:
        """
        python3 /beegfs/home/bbein/polymerase_snakemake_wrapper/prepare_kmc.py \
        -i {input.kmc} \
        -o {output.script};
        kmc_tools \
        complex \
        {output.script} \
        2> {log}
        """

## Get differences for all polymerases
rule kmc_diff:
    input:
        kmc = expand("kmc/{sample}.res.kmc_suf", sample = samples),
    output:
        script = "kmc/{sample}_diff_script",
        kmc_out_pre = "kmc/{sample}_diff.kmc_pre",
        kmc_out_suf = "kmc/{sample}_diff.kmc_suf"
    threads:
        10
    resources:
        mem_mb  = 40000
    params:
        sample = "{sample}",
        memory=40
    log:
        "logs/kmc_diff/{sample}.log",
    envmodules:
        "kmc/3.2.4"
    shell:
        """
        python3 /beegfs/home/bbein/polymerase_snakemake_wrapper/prepare_kmc.py \
        -i {input.kmc} \
        -r {params.sample};
        kmc_tools \
        complex \
        {output.script} \
        2> {log}
        """

# rule kmc_stats: