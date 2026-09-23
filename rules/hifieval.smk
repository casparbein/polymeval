minimap_wrapper = f"{wrapper_versions['minimap']}/bio/minimap2/aligner"
seqtk_wrapper = f"{wrapper_versions['seqtk']}/bio/seqtk"

CORRECTED = {
    "hifiasm": "assemblies/hifiasm/{sample}.ec.fa",
    "lja":     "assemblies/lja/{sample}/k5001/corrected_reads.fasta",
    "verkko":  "assemblies/verkko/{sample}/hifi-corrected.fasta.gz",
}

def corrected_reads(wc):
    sample, corrector = (wc.asm_id.rsplit("__", 1) if "__" in wc.asm_id
                         else (wc.asm_id, ASSEMBLERS[0]))
    return CORRECTED[corrector].format(sample=sample)

## Whether error correction includes HP compression
HPC_CORRECTOR = {"hifiasm": False, "lja": True, "verkko": True}

def corrector_of(asm_id):
    return asm_id.rsplit("__", 1)[1] if "__" in asm_id else ASSEMBLERS[0]

def hifieval_target(wc):
    return (f"hifieval/hpc/{wc.asm_id}.target.fa"
            if HPC_CORRECTOR[corrector_of(wc.asm_id)]
            else f"assemblies/{wc.asm_id}.fa")

def hifieval_raw_query(wc):
    s = asm_sample(wc.asm_id)
    return (f"hifieval/hpc/{s}.raw.fa"
            if HPC_CORRECTOR[corrector_of(wc.asm_id)]
            else f"raw_reads/{s}.fastq.gz")

rule hpc_target:
    input:  
        "assemblies/{asm_id}.fa"
    output: 
        temp("hifieval/hpc/{asm_id}.target.fa")
    conda:  
        "../envs/seqtk.yaml"
    log:    
        "logs/hpc_target/{asm_id}.log"
    params:
        command="hpc",
    wrapper:  
        seqtk_wrapper

rule hpc_raw_reads:
    input:  
        "raw_reads/{sample}.fastq.gz"
    output: 
        temp("hifieval/hpc/{sample}.raw.fa")
    conda:  
        "../envs/seqtk.yaml"
    log:    
        "logs/hpc_raw_reads/{sample}.log"
    params:
        command="hpc",
    wrapper:
        seqtk_wrapper

## Approximation of read error stats with hifieval
## Align raw reads
rule hifieval_align_raw:
    input:
        target=hifieval_target,
        query=hifieval_raw_query,
    output:
        temp("alignments/{asm_id}.raw.paf"),
    log:
        "logs/hifieval_align_raw/{asm_id}.log",
    params:
        extra="-cx map-hifi --secondary=no --paf-no-hit --cs", 
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 20
    resources:
        mem_mb = 100000
    wrapper:
       minimap_wrapper

## Aligned error corrected reads
rule hifieval_align_ec:
    input:
        target=hifieval_target,
        query=corrected_reads,
    output:
        temp("alignments/{asm_id}.ec.paf"),
    log:
        "logs/hifieval_align_ec/{asm_id}.log",
    params:
        extra="-cx map-hifi --secondary=no --paf-no-hit --cs", 
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 20
    resources:
        mem_mb = 100000
    wrapper:
        minimap_wrapper

## Normally optional: If read file contains empty reads (Revio demo from PacBio does), hifieval will choke on this:
rule remove_empty:
    input:
        ec = "alignments/{asm_id}.ec.paf",
        raw = "alignments/{asm_id}.raw.paf",
    output:
        ec = temp("alignments/{asm_id}.ec.clean.paf"),
        raw = temp("alignments/{asm_id}.raw.clean.paf"),
    log: 
        "logs/remove_empty/{asm_id}.log",
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        awk '$2 > 0 {{print}}' {input.ec}  > {output.ec};
        awk '$2 > 0 {{print}}' {input.raw} > {output.raw};
        """

## Run Hifieval
rule hifieval_compare:
    input:
        raw ="alignments/{asm_id}.raw.clean.paf", 
        ec ="alignments/{asm_id}.ec.clean.paf", 
    output:
        metric = "hifieval/{asm_id}.metric.eval.tsv",
        rdl_eval = "hifieval/{asm_id}.rdlvl.eval.tsv",
        summary = temp("hifieval/{asm_id}.summary.tsv"),
    params:
        out_base = "hifieval/{asm_id}",
    threads:
        1
    resources:
        mem_mb = 100000,
    log:
        "logs/hifieval_compare/{asm_id}.log",
    conda:
        "../envs/hifieval.yaml",
    shell:
        """
        hifieval.py \
        -o {params.out_base} \
        -r {input.raw} \
        -c {input.ec} \
        2> {log}
        """