## Genomescope 
rule genomescope:
    input:
        hist="out/stats/{sample}.kmc_hist.txt",
    output:
        multiext(
            "genomescope/{sample}/",
            "linear_plot.png",
            "log_plot.png",
            "model.txt",
            "progress.txt",
            "SIMULATED_testing.tsv",
            "summary.txt",
            "transformed_linear_plot.png",
            "transformed_log_plot.png",
        ),
    log:
        "logs/genomescope/{sample}.log",
    params:
        extra="--kmer_length 25 --testing -m 50000",
    wrapper:
        f"{wrapper_versions['genomescope']}/bio/genomescope"