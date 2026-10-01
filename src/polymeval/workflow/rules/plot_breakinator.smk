## Plot breakinator output
rule plot_breakinator:
    input:
        expand("out/files/{sample}_breakinator_summary.txt", sample = samples),
    output:
        plot = "out/plots/breakinator.pdf",
    threads:
        1
    resources:
        mem_mb = 50000
    params:
        chimera_path = "out/files/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples),
    conda:
        "../envs/plot_chimeras.yaml"
    log:
        "logs/plot_breakinator/plot_breakinator.log",
    script:
        "../scripts/plot_breakinator.R"