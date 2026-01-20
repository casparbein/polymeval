## Plot GC vs Norm Depth for all polymerases
rule plot_gc:
    input:
        pandepth = expand("out/files/{sample}_pandepth_all.win.stat.gz", sample = samples),
    output:
        boxplot = "out/plots/GC_boxes.pdf",
        mean_line = "out/plots/GC_mean.pdf",
        median_line = "out/plots/GC_median.pdf",
        iqr_line = "out/plots/GC_iqr.pdf",
        lci_bar = "out/plots/LCI_bar.pdf",
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
        "logs/plot_gc/plot_gc.log"
    script:
        "../scripts/gc_bins.R"