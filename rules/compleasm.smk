rule run_compleasm:
    input:
        "assemblies/{asm_id}.fa"
    output:
        "compleasm/{asm_id}_compleasm/summary.txt",
        temp(directory(f"compleasm/{{asm_id}}_compleasm/{config['compleasm_db']}/hmmer_output")),
    threads:
        10
    resources:
        mem_mb=50000
    params:
        outname = "compleasm/{asm_id}_compleasm",
        database = config["compleasm_db"],
        database_path = config["compleasm_db_path"],
    log:
        "logs/run_compleasm/{asm_id}.log"
    conda:
        "../envs/compleasm.yaml"
    shell:
        """
        compleasm \
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
        summary = "compleasm/{asm_id}_compleasm/summary.txt",
    output:
        "compleasm/{asm_id}_summary.rf.txt"
    log:
        "logs/reformat_compleasm/{asm_id}.log"
    shell:
        """
        cat {input.summary} | sed -e 's/:/\t/g' -e 's/%, /\t/g' | head -n6 | tail -n5 > {output} 2> {log}
        """