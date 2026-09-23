LIB     = config["compleasm_db_path"]
LINEAGE = config["compleasm_db"]
HAVE_LINEAGE = os.path.isdir(os.path.join(LIB, LINEAGE))

localrules: compleasm_download

## Download compleasm lib if it is not available
rule compleasm_download:
    output: 
        touch(os.path.join(LIB, ".polymeval." + LINEAGE + ".ok"))
    params:
        lib     = LIB,
        lineage = LINEAGE.rsplit("_odb", 1)[0],
        odb     = "odb" + LINEAGE.rsplit("_odb", 1)[1],
    retries: 3
    log: 
        "logs/compleasm_download/download.log"
    conda:
        "../envs/compleasm.yaml"
    shell:
        """
        compleasm download \
        -L {params.lib} \
        --odb {params.odb} \
        {params.lineage} &> {log}
        """

## Run compleasm
rule run_compleasm:
    input:
        "assemblies/{asm_id}.fa",
        lib = [] if HAVE_LINEAGE else os.path.join(LIB, ".polymeval." + LINEAGE + ".ok"),
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
        -L {input.lib} \
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