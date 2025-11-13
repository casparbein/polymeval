rule run_compleasm:
    input:
        "assemblies/{sample}.fa"
    output:
        "compleasm/{sample}_compleasm/summary.txt",
        temp(directory(f"compleasm/{{sample}}_compleasm/{config['compleasm_db']}/hmmer_output")),
    threads:
        10
    resources:
        mem_mb=50000
    params:
        outname = "compleasm/{sample}_compleasm",
        database = config["compleasm_db"],
        database_path = config["compleasm_db_path"],
    log:
        "logs/run_compleasm/{sample}.log"
    conda:
        "../envs/compleasm.yaml"
    envmodules:
        "compleasm/0.2.7"
    shell:
        """
        compleasm.py \
        run \
        -a {input} \
        -o {params.outname} \
        -l {params.database} \
        -L {params.database_path} \
        -t {threads} \
        2> {log}
        """

## reformat stats so they can be read in easily in R
rule reformat_compleasm:
    input:
        "compleasm/{sample}_compleasm/summary.txt"
    output:
        "compleasm/{sample}_summary.rf.txt"
    log:
        "logs/reformat_compleasm/{sample}.log"
    shell:
        """
        cat {input} | sed -e 's/:/\t/g' -e 's/%, /\t/g' | head -n6 | tail -n5 > {output} 2> {log}
        """