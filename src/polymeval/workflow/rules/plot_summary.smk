## Plot GC vs Norm Depth for all polymerases
rule plot_summary:
    input:
        pandepth = expand("out/files/{sample}_pandepth_all.win.stat.gz", sample = samples),
        pafs = expand("alignments/{sample}_to_ref.paf", sample = samples),
        breakpoint_bed = "out/files/breakpoints.bed",
    output:
        summary_plot = "out/plots/summary_plot.pdf",
    threads:
        1
    resources:    
        mem_mb = 100000
    params:
        pandepth_path = "out/files/",
        paf_path = "alignments/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples),
    conda:
        "../envs/plot_summary.yaml"
    log:
        "logs/plot_summary/plot_summary.log"
    script:
        "../scripts/summary_view.R"