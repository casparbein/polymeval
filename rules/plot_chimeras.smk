## Extract chimeras from reads
rule plot_chimeras:
    input:
        map_type = expand("out/files/{sample}_mappingType.txt", sample = samples),
        prim_hist = expand("out/files/{sample}_prim.hist.txt", sample = samples),
        supp_hist = expand("out/files/{sample}_suppl.hist.txt", sample = samples),
        rl_hist = expand("out/files/{sample}_rl.hist.txt", sample = samples),
    output:
        plot = "out/plots/chimeras.pdf",
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
        "logs/plot_chimeras/plot_chimeras.log",
    script:
        "../scripts/plot_chimeras.R"