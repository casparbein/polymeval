## Create final output plots and stats
rule final_output_stats:
    input:
        faidx_files = expand("assemblies/{sample}.fa.fai", sample = samples),
        compleasm_files = expand("compleasm/{sample}_summary.rf.txt", sample = samples),
        merqury_files = expand("merqury/{sample}_slf/{sample}_slf.qv", sample = samples),
        hifieval_files = expand("hifieval/{sample}.summary.tsv", sample = samples) if config["hifieval"] else [],
        seqkit = "out/stats/seqkit_all.tsv"
    output:
        ng_table = "out/stats/all.NG.table.txt",
        compleasm_table = "out/stats/all.compleasm.stats.txt",
        final_plot = "out/stats/all.summary_plot.pdf",
        hifieval_table = "out/stats/hifieval_table.tsv" if config["hifieval"] else [],
        merqury_table = "out/stats/merqury_table.tsv",
    params:
        faidx_path = "assemblies/",
        compleasm_database = config["compleasm_db"],
        compleasm_path = "compleasm/",
        hifieval_path = "hifieval/" if config["hifieval"] else [],
        merqury_path = "merqury/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples),
    threads:
        1
    resources:
      mem_mb = 100000
    log:
        "logs/final_output_stats/final.log",
    conda:
        "../envs/ref_free_summary.yaml"
    script:
        "../scripts/ref_free_summary.R"

## Create read output plots (?) and stats
rule read_output_stats:
    input:
        readeval_dump = expand("out/stats/{sample}.rdeval_dump.tsv",sample = samples),
        qchist = expand("out/stats/{sample}.qchist.txt", sample = samples),
    output:
        qv_plot = "out/stats/qv_plot.pdf",
        standard_hist_plot_readeval = "out/stats/standard_hist_plot_readeval.pdf",
        gc_hist_hist = "out/stats/gc_hist_hist.pdf",
        rdeval_like = "out/stats/rdeval_like.pdf",
    params:
        stats_path = "out/stats/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples, allow_missing=True),
    threads:
        1
    resources:
        mem_mb = 50000
    log:
        "logs/read_output_stats/final.log",
    conda:
        "../envs/ref_free_summary.yaml"
    script:
        "../scripts/read_stats.R"