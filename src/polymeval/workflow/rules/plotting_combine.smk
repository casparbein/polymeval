## Create final output plots and stats
rule final_output_stats_combine:
    input:
        faidx_files = expand("assemblies/{out}.fa.fai", out = out_assembly_names),
        compleasm_files = expand("compleasm/{out}_summary.rf.txt", out = out_assembly_names),
        merqury_files = expand("merqury/{out}_slf/{out}_slf.qv", out = out_assembly_names),
        hifieval_files = expand("hifieval/{out}.summary.tsv", out = out_assembly_names) if config["hifieval"] else [],
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
        sample_names = expand("{out}", out = out_assembly_names),
    threads:
        1
    resources:
      mem_mb = 100000
    log:
        "logs/final_output_stats/final.log",
    conda:
        "../envs/ref_free_summary.yaml"
    script:
        "../scripts/ref_free_summary_combine.R"