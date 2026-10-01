## Get chromosomal breakpoints against best reference
rule get_breakpoints:
    input:
        pafs=expand("alignments/{sample}_to_ref.paf", sample = samples),
    output:
        bed = "out/files/breakpoints.bed",
        plot = "out/plots/breakpoint_matrix.pdf"
    threads:
        1
    resources:
        mem_mb = 100000
    params:
        paf_path = "alignments/",
        colors = config["colors"],
        sample_names = expand("{sample}", sample = samples),
    conda:
        "../envs/get_breakpoints.yaml"
    log:
        "logs/get_breakpoints/breakpoints.log",
    script:
        "../scripts/calculate_contig_breaks.R"