## map single_polymerase (or any other) assemblies
rule minimap2_paf:
    input:
        target=reference_seq,
        query="assemblies/{sample}.fa",
    output:
        "alignments/{sample}_to_ref.paf",
    log:
        "logs/minimap2_paf/{sample}.log",
    params:
        extra="-x asm5 -c -eqx -secondary=no",  # optional
        sorting="",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    resources:
        mem_mb = 100000
    threads: 20
    wrapper:
        "v7.6.0/bio/minimap2/aligner"