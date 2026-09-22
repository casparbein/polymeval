samtools_wrapper = f"{wrapper_versions['samtools']}/bio/samtools/faidx"
hifiasm_wrapper = f"{wrapper_versions['hifiasm']}/bio/hifiasm"

## Where each assembler leaves its final FASTA. This table is the only place those
## paths are recorded; adding an assembler means adding one line here and one rule below.
ASSEMBLY_FASTA = {
    "hifiasm": "assemblies/hifiasm/{sample}.fa",
    "flye":    "assemblies/flye/{sample}/assembly.fasta",
    "lja":     "assemblies/lja/{sample}/assembly.fasta",
    "verkko":  "assemblies/verkko/{sample}/assembly.fasta",
}

def assembly_source(wc):
    sample, assembler = (wc.asm_id.rsplit("__", 1) if "__" in wc.asm_id
                         else (wc.asm_id, ASSEMBLERS[0]))
    return ASSEMBLY_FASTA[assembler].format(sample=sample)

## Polymerase-specific assembly with hifiasm
rule hifiasm:
    input:
        fasta="raw_reads/{sample}.fastq.gz",
    output:
        "assemblies/hifiasm/{sample}.p_ctg.gfa",
        temp("assemblies/hifiasm/{sample}.ec.fa") if config["hifieval"] else [],
        temp("assemblies/hifiasm/{sample}.r_utg.gfa"),
        temp("assemblies/hifiasm/{sample}.p_utg.gfa"),
        temp("assemblies/hifiasm/{sample}.a_ctg.gfa"),
        temp("assemblies/hifiasm/{sample}.a_ctg.lowQ.bed"),
        temp("assemblies/hifiasm/{sample}.a_ctg.noseq.gfa"),
        temp("assemblies/hifiasm/{sample}.p_ctg.lowQ.bed"),
        temp("assemblies/hifiasm/{sample}.p_ctg.noseq.gfa"),
        temp("assemblies/hifiasm/{sample}.p_utg.lowQ.bed"),
        temp("assemblies/hifiasm/{sample}.p_utg.noseq.gfa"),
        temp("assemblies/hifiasm/{sample}.r_utg.lowQ.bed"),
        temp("assemblies/hifiasm/{sample}.r_utg.noseq.gfa"),
        temp("assemblies/hifiasm/{sample}.ec.bin"),
        temp("assemblies/hifiasm/{sample}.ovlp.reverse.bin"),
        temp("assemblies/hifiasm/{sample}.ovlp.source.bin"),
    log:
        "logs/hifiasm/{sample}.log",
    params:
        extra=f"--primary -l 3 --write-ec --hg-size {config['hg_size']}" if config["hifieval"] and config["hg_size"] else " --primary -l 3  --write-ec " if config["hifieval"] and not config["hg_size"] else  f" --primary -l 3 --hg-size {config['hg_size']}" if not config["hifieval"] and config["hg_size"] else "--primary -l 3",
    threads: 32
    resources:
        mem_mb=200000,
    wrapper:
       hifiasm_wrapper

rule get_fasta:
    input:  "assemblies/hifiasm/{sample}.p_ctg.gfa"
    output: "assemblies/hifiasm/{sample}.fa"
    log:    "logs/get_fasta/{sample}.log"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} > {output} 2> {log}
        """

## Flye assembly
rule flye:
    input:  "raw_reads/{sample}.fastq.gz"
    output: fa = "assemblies/flye/{sample}/assembly.fasta",
            gfa = "assemblies/flye/{sample}/assembly_graph.gfa",
    params:
        outdir = "assemblies/flye/{sample}",
        gsize  = f"--genome-size {config['hg_size']}" if config.get("hg_size") else "",
        extra  = config.get("flye_extra", ""),
    threads: 32
    resources: mem_mb = 200000
    conda: "../envs/flye.yaml"
    log: "logs/flye/{sample}.log"
    shell:
        """
        flye \
        --pacbio-hifi {input} \
        --out-dir {params.outdir} \
        --threads {threads} \
        {params.gsize} \
        {params.extra} &> {log}
        """

## LJA assembly
rule lja:
    input:  
        "raw_reads/{sample}.fastq.gz"
    output: 
        "assemblies/lja/{sample}/assembly.fasta",
        temp("assemblies/lja/{sample}/01_TopologyBasedCorrection/corrected_reads.fasta")
    params:
        outdir  = "assemblies/lja/{sample}",
        diploid = "--diploid" if config.get("lja_diploid") else "",
    threads: 32
    resources: mem_mb = 200000
    conda: "../envs/lja.yaml"
    log: "logs/lja/{sample}.log"
    shell:
        """
        lja \
        -o {params.outdir} \
        --reads {input} \
        -t {threads} \
        --diploid &> {log}
        """

## Verkko assembly (might not work since it is local)
rule verkko:
    input:  
        "raw_reads/{sample}.fastq.gz"
    output: 
        "assemblies/verkko/{sample}/assembly.fasta",
        temp("assemblies/verkko/{sample}/hifi-corrected.fasta.gz")
    params:
        outdir  = "assemblies/verkko/{sample}",
        mem_gb  = lambda wc, resources: max(1, resources.mem_mb // 1024),
        extra   = config.get("verkko_extra", ""),
    threads: 32
    resources: mem_mb = 100000
    conda: "../envs/verkko.yaml"
    log: "logs/verkko/{sample}.log"
    shell:
        """
        verkko \
        -d {params.outdir} \
        --hifi {input} \
        --local \
        --local-cpus {threads} \
        --local-memory {params.mem_gb} \
        {params.extra} &> {log}
        """

rule collect_assembly:
    input:  assembly_source
    output: "assemblies/{asm_id}.fa"
    run:
        os.symlink(os.path.relpath(input[0], os.path.dirname(output[0])), output[0])

rule samtools_faidx:
    input:  
        "assemblies/{asm_id}.fa",
    output: 
        "assemblies/{asm_id}.fa.fai",
    log:    
        "logs/faidx/{asm_id}.log",
    params:
        extra="",
    wrapper:
        samtools_wrapper