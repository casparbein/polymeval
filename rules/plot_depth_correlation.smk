## Plot depth correlations between all polymerases
rule plot_depth_correlation:
    input:
        pandepth = expand("out/files/{sample}_pandepth_all.win.stat.gz", sample = samples),
        contig_breaks = "out/files/breakpoints.bed"
    output:
        overall_matrix = "out/plots/overall_depth_matrix.pdf",
        breakpoint_matrix = "out/plots/breakpoint_depth_matrix.pdf",
        by_input = "out/plots/dropout_summary.pdf"
    threads:
        1
    resources:    
        mem_mb = 100000
    params:
        pandepth_path = "out/files/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples),
    conda:
        "../envs/plot_depth_correlation.yaml"
    log:
        "logs/plot_depth_correlation/plot_depth_correlation.log"
    script:
        "../scripts/depth_correlation.R"