import argparse
import logging
import sys
import os
from ruamel.yaml import load_all, YAML, comments
from ruamel.yaml.scalarstring import SingleQuotedScalarString, DoubleQuotedScalarString
import subprocess
import shutil
import signal

## import helper script
import get_downsample_rates

## Logging
logger = logging.getLogger("polymeval")
logging.basicConfig(level=logging.INFO, format="[%(levelname)-8s] %(message)s")

## Absolute paths on file system
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "."))

## Suffixes of input read files:
CANONICAL_SUFFIX = {
    ".fq":        ".fastq",
    ".fq.gz":     ".fastq.gz",
    ".dup.fq":    ".dup.fastq",
    ".dup.fq.gz": ".dup.fastq.gz",
    ".fasta":     ".fa",
    ".fna":       ".fa",
}

ASM_SUFFIXES = (".fa", ".fasta", ".fna")

## Helper functions
def split_suffix(name, suffixes):
    """Return (basename, matched_suffix), or (None, None) if nothing matches."""
    for suf in suffixes:
        if name.endswith(suf):
            return name[: -len(suf)], suf
    return None, None

## Helper functions for mounting apptainer paths:
def get_mount_point(path):
    path = os.path.realpath(os.path.abspath(path))
    while not os.path.ismount(path):
        path = os.path.dirname(path)
    return path

def get_apptainer_bind_args(paths):
    mount_points = set()
    for path in paths:
        if path and os.path.exists(path):
            mp = get_mount_point(path)
            if mp != '/':  # skip root
                mount_points.add(mp)
    return " ".join(f"-B {mp}:{mp}" for mp in mount_points)


def link_inputs(src_path, dst_path, suffixes, down_list=None,
                skip_bases_ending=None, kind="input"):
    """Link every accepted file from src_path into dst_path under its canonical name.

    suffixes           accepted source extensions; matched longest-first
    down_list          restrict to these sample basenames, or None for all
    skip_bases_ending  ignore files whose basename ends with this (".ec", ".dup")
    kind               noun used in error messages

    Returns [(basename, canonical_suffix), ...] sorted by basename.
    """
    os.makedirs(dst_path, exist_ok=True)
    real_src = os.path.realpath(os.path.abspath(src_path))
    suffixes = tuple(sorted(suffixes, key=len, reverse=True))

    found = {}
    for name in sorted(os.listdir(real_src)):
        base, suf = split_suffix(name, suffixes)
        if base is None:
            continue
        if skip_bases_ending and base.endswith(skip_bases_ending):
            continue
        if down_list is not None and base not in down_list:
            continue
        if base in found:
            logger.critical("Two %s files map to sample %r in %s: %s and %s. Keep one.",
                            kind, base, src_path, found[base][1], name)
            sys.exit(1)
        found[base] = (CANONICAL_SUFFIX.get(suf, suf), name)

    # nothing is linked until the whole directory has been accepted
    for base, (canonical, name) in found.items():
        dst_file = os.path.join(dst_path, base + canonical)
        if not os.path.lexists(dst_file):
            os.symlink(os.path.join(real_src, name), dst_file)

    return sorted((b, c) for b, (c, _) in found.items())

## Check whether compleasm lib is available, otherwise, write to user cache
def resolve_compleasm_lib(arg):
    if arg:
        if not os.access(os.path.expanduser(arg), os.W_OK):
            logger.critical("compleasm needs write access to %s (it refreshes file_versions.tsv on every "
                    "run, even when the lineage is already present). Copy the lineage to a "
                    "writable directory and pass it with --compleasm_db_path, set "
                    "POLYMEVAL_COMPLEASM_LIBS, or do not set this parameter, leading to automatic download "
                    "of the specified db to the user's cache. ", os.path.expanduser(arg))
            sys.exit(1)
        else:
            return os.path.abspath(os.path.expanduser(arg))
    env = os.environ.get("POLYMEVAL_COMPLEASM_LIBS")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    cache = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(cache, "polymeval", "compleasm")


def link_reads(src, dst, suffixes, down_list=None, skip_bases_ending=None):
    return link_inputs(src, dst, suffixes, down_list, skip_bases_ending, kind="read")

def link_assemblies(src, dst, down_list=None):
    return link_inputs(src, dst, ASM_SUFFIXES, down_list,
                       skip_bases_ending=".ec", kind="assembly")

def get_snakefile_path(name="Snakefile"):
    snakefile = os.path.join(base_dir, name)
    return snakefile

def get_cluster_configfile_path(name="config.yaml"):
    cluster_configfile =  os.path.join(base_dir, "prof")
    return cluster_configfile

def get_kmc_prep_path(name="prepare_kmc.py"):
    kmc_prep_file =  os.path.join(base_dir, "scripts", name)
    return kmc_prep_file

def format_list(in_list):
    format_list = comments.CommentedSeq(in_list)
    format_list.fa.set_flow_style()
    return format_list

def run_snakemake(snake_file, 
                  config_file, 
                  conda_path,
                  directory_name="polymeval_test",
                  dryrun = True, 
                  snake_default = False, 
                  rerun_triggers = False, 
                  working_dir=os.getcwd(),
                  updated_rule = False,
                  use_conda = True,
                  local_run = False,
                  use_apptainer = False,
                  apptainer_args = None):
    
    ## move to the working directory
    current_wd =  os.getcwd()
    if working_dir == current_wd:
        snakemake_dir = os.path.join(current_wd, directory_name)
    else:
        snakemake_dir = os.path.join(working_dir, directory_name)

    if not os.path.exists(snakemake_dir):
        os.makedirs(snakemake_dir)
    
    ## Change to snakemake dir as wd
    os.chdir(snakemake_dir)
    
    ## Run snakemake from the command line
    cmd = ['snakemake']

    snake_file = get_snakefile_path(snake_file)
    cmd += ['--snakefile', snake_file]
    
    if not local_run:
        config_file = get_cluster_configfile_path("config.yaml")
        cmd += ['--profile', config_file]
    else:
        ## Add these arguments to CLI so users can set them
        cmd += ['--cores', '4', '--max-threads', '4', '--resources', 'mem_mb=30000']

    cmd += ['--directory', snakemake_dir]

    if use_conda and not conda_path:
        conda_path = os.path.join(base_dir, ".snakemake/conda")
        cmd += ['--use-conda', '--conda-prefix', conda_path]
    elif use_conda and conda_path:
        cmd += ['--use-conda', '--conda-prefix', conda_path]
    elif not use_conda and conda_path:
        logger.critical("--use-conda was not specified, but a reference conda path given. Since use_conda is activated by default, do not activate it if you also pass a conda path")
        sys.exit(1)
        #sys.exit("--use-conda was not specified, but a reference conda path given. Since use_conda is activated by default, do not activate it if you also pass a conda path")
    elif not use_conda and not conda_path:
        logger.warning("WARNING: --use-conda is deactivated and no conda path is found. Most likely the pipeline will fail, unless all tools are installed on the user's machine and availabe in $PATH")

    if use_apptainer:
        apptainer_path = os.path.join(base_dir, ".snakemake/singularity")
        cmd += ['--use-apptainer', '--apptainer-prefix', apptainer_path]
        if apptainer_args:
            cmd += ['--apptainer-args', apptainer_args]

    #if use_apptainer:
    #    apptainer_path = os.path.join(base_dir, ".snakemake/singularity")
    #    cmd += ['--use-apptainer','--apptainer-prefix', apptainer_path, '--apptainer-args', f'"{apptainer_args}"']

    if dryrun:
        cmd.append('--dry-run')
    
    if rerun_triggers:
        cmd += ['--rerun-triggers', 'mtime']    
    
    if snake_default:
        default_snakemake_args = ["--rerun-incomplete", "--keep-going"]
        cmd += default_snakemake_args

    #print_cmd = " ".join(cmd)
    #print("The following command will be run for the polymeval pipeline: {}".format(print_cmd))

    ## Run pipeline proper
    #try:
    #    result1 = subprocess.run(cmd, check=True)
    #    print(result1.stdout)
    #except subprocess.CalledProcessError as e:
    #    print(f"Command failed with return code {e.returncode}")

    print_cmd = " ".join(cmd)
    logger.info("The following command will be run for the polymeval pipeline: %s", print_cmd)


    result1 = subprocess.Popen(cmd,preexec_fn=os.setpgrp)
    try:
        # Wait for Snakemake to finish normally
        result1.wait()
    except KeyboardInterrupt:
        # 2. Wrapper catches the Ctrl+C
        logger.info("Intercepted Ctrl+C. Forwarding to Snakemake...")
        
        # 3. Send SIGINT specifically to the Snakemake process group
        # We use the negative PID to signal the entire group (Snakemake + its workers)
        os.killpg(result1.pid, signal.SIGINT)
        
        logger.info("Waiting for Snakemake to finish Slurm/cleanup (up to several minutes)...")
        
        # 4. Stay here until Snakemake exits. 
        # We ignore further KeyboardInterrupts during this specific wait.
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        result1.wait()
        
        logger.info(f"Snakemake has exited (Code: {result1.returncode}).")
    sys.exit(result1.returncode)


DESCRIPTION = '''
polymeval - a snakemake pipeline to streamline benchmarking and comparing PacBio HiFi datasets amplified with different polymerases.
'''

def argument_parser():
    """Parse CMD args."""
    app = argparse.ArgumentParser(description=DESCRIPTION)
    #app = argparse.ArgumentParser(description=DESCRIPTION)

    ## In which mode should the pipeline be run
    run_mode = app.add_mutually_exclusive_group(
                                            required=True
                                            )
    run_mode.add_argument(
    "-s",
    "--standard", 
    action="store_true", 
    dest="standard",
    help=
    """Run the standard pipeline (read stats, hifiasm assembly, compleasm stats, merqury on a set of input reads).
    """
    )

    run_mode.add_argument(
    "-d",
    "--downsample", 
    action="store_true",
    dest="downsample",
    help=
    """Run the downsample pipeline (Must provide a file to a seqkit stats output. 
    By default runs downsampling, hifiasm assembly, compleasm stats, merqury)
    """
    )

    run_mode.add_argument(
    "-c",
    "--combine", 
    action="store_true",
    dest="combine",
    help=
    """combine/downsample different read sets to see whether they complement each other.
    By default takes a set of up to 5 read sets, finds the smallest, downsamples to 1/n of that number of nucleotides
    and creates hifiasm assembly, compleasm stats, merqury for those assemblies.
    """
    )

    run_mode.add_argument(
    "-r",
    "--reference", 
    action="store_true",
    dest="reference",
    help=
    """Evaluate read sets based on a reference assembly. Per-read assemblies also have to be provided, those can be created by 
    any of the other polymeval modes. This mode will evaluate raw reads and per-polymerase assemblies based on a reference, and
    is therefore supposed to be executed downstream of standard.
    """
    )

    run_mode.add_argument(
    "-v",
    "--variant_calling_benchmarks", 
    action="store_true",
    dest="variant_calling",
    help=
    """IMPORTANT: Only works for human data sequenced from the HG002 sample.
    Call variants of (amplified) HiFi reads against the hg38 (hg37 for structural variants) with deepvariant (sniffles) and
    evaluate against benchmark sets with happy (truvari). If structural variant accuracy should also be evaluated, set 
    --structural_variant_calling flag. All reference benchmark files must be in a path that can be set with --benchmark_path.
    """
    )

    run_mode.add_argument(
    "-fb",
    "--fetch_benchmarks", 
    action="store_true",
    help=
    '''Download and prepare all human variant benchmark data into --benchmark_path, then exit.
    This can be run before a variant benchmarking run is started to have all necessary ground truth files present and in the right format.
    ''')

    app.add_argument(
    "-A", 
    "--assembler",
    action="store", 
    dest="assemblers", 
    default="hifiasm",
    help=
    """Comma-separated list of assemblers: hifiasm, flye, lja, verkko.
    Default: hifiasm. Multiple assemblers are run on every sample and
    compared in the summary plot.
    """)

    app.add_argument(
    "-pw",
    "--pairwise", 
    action="store_true",
    dest="pairwise",
    help=
    """For combine mode: If only pairwise combinations should be run (instead of the default: all combos). Allows up to 8 input samples
    instead of 5 in default combine.
    """
    )

    app.add_argument(
    "-i", 
    "--input_reads",
    action="store",
    dest="in_reads",
    type=str,
    help=
    '''Path to directory with input reads. It is strongly recommended to give the 
    reads meaningful names (those names will be propagated down to the output files). 
    Also, including '.' or '-' characters in the name is not supported 
    (Example: sample1.1.fastq.gz or sample1-1.fastq.gz does not work, but sample1.fastq.gz does).
    ''')

    app.add_argument(
    "-du", 
    "--remove_dups",
    action="store_true",
    dest="remove_dups",
    help=
    '''Whether PacBio's pbmarkdup should be run (Here, PCR dups are automatically marked and removed).
    Input reads in the path passed to -i/--input_reads have to be named SAMPLE.dup.fastq.gz instead of 
    just SAMPLE.fastq.gz. Only works in standard mode.
    ''')

    app.add_argument(
    "-ia", 
    "--input_assemblies",
    action="store",
    dest="in_assemblies",
    type=str,
    help=
    '''For 'reference' mode: Path to directory with input assemblies. Naming as with input reads, and names have to match.
    Normally, these assemblies will have been created with --standard, --combine or --downsample and the naming therefore
    will work out naturally.
    ''')

    app.add_argument(
    "-rsq", 
    "--reference_sequence",
    action="store",
    dest="reference_seq",
    default="",
    type=str,
    help=
    '''Path to reference sequence (Normally a reference quality assembly) to which the different datasets should be compared,
    in fasta format.
    ''')

    app.add_argument(
    "-db", 
    "--compleasm_db",
    action="store",
    dest="compleasm_db",
    type=str,
    help=
    '''Name of the compleasm db that should be used to compute assembly completness. 
    Must be present in the specified path, otherwise will be downloaded automatically,
    see --compleasm_db_path
    ''')

    app.add_argument(
    "-dbp", 
    "--compleasm_db_path",
    action="store",
    dest="compleasm_db_path",
    default=None,
    type=str,
    help=
    '''Path to accessible odb library stored centrally on the HPC. If not provided, 
    compleasm will download libraries on the fly and put it on the user's cache, which takes time and disk space.
    ''')

    app.add_argument(
    "-co", 
    "--colors",
    action="store",
    dest="colors",
    default="",
    type=str,
    help=
    '''Path to a tsv file listing sample names in column 1 and 
    desired colors (in hexadecimal code or R notation) in column2. Can be used to add colors to the output plots.
    If not provided, will automatically create a color scale in R.
    ''')

    app.add_argument(
    "-hi", 
    "--hifieval",
    action="store_true",
    dest="hifieval",
    help=
    '''Turn on hifieval
    ''')

    app.add_argument(
    "-hg", 
    "--hg_size",
    action="store",
    dest="hg_size",
    help=
    '''Estimated size for homozygous genome length, for hifiasm's --hg-size.
    Should be in k,m or g (for instance: 2g or 700m).
    ''')

    app.add_argument(   
    "--lja_path", 
    action="store", 
    dest="lja_path", 
    default=None,
    help=
    '''Absolute path to an LJA binary.
    ''')

    ## Define exactly what the output files would be here
    app.add_argument(
    "-rd", 
    "--readstats",
    action="store_true",
    dest="readstats",
    help=
    '''Turn on read stats computation:
    Read length histogram, QC histogram, rdeval dump and associated plots
    ''')

    app.add_argument(
    "-k", 
    "--kmc",
    action="store_true",
    dest="kmc",
    help=
    '''Turn on KMC kmer counting. Genomescope will be run on resulting kmer histograms.
    ''')

    # app.add_argument(
    # "-pa", 
    # "--pandepth",
    # action="store_true",
    # dest="pandepth",
    # help=
    # '''Turn on Pandepth depth mapping. Most reference mode analyses will need this.
    # ''')

    # app.add_argument(
    # "-pp", 
    # "--pandepth_path",
    # action="store",
    # dest="pandepth_path",
    # default = "",
    # help=
    # '''In case pandepth is installed but not in the user's $PATH, provide absolute path to pandepth.
    # ''')

    app.add_argument(
    "-km", 
    "--kmer_length",
    action="store",
    dest="kmer_length",
    default=25,
    help=
    '''K-mer length used by meryl and KMC (if enabled). Default is 25.
    ''')

    app.add_argument(
    "-se", 
    "--seed",
    action="store",
    dest="rasusa_seed",
    default=100,
    help=
    '''Seed for rasusa downsampling.
    ''')

    app.add_argument(
    "-bp", 
    "--benchmark_path",
    action="store",
    dest="benchmark_path",
    default=None,
    help=
    '''Path to where benchmark files for human variant calling are stored.
    ''')
    
    app.add_argument(
    "--benchmark_releases", 
    default="v5.0q,cmrg,tandem_repeats",
    help=
    '''Comma-separated GIAB releases to fetch: v5.0q, v4.2.1, cmrg, tandem_repeats, NIST_SV_v0.6.
    ''')

    app.add_argument(
    "-svc", 
    "--structural_variant_calling",
    action="store_true",
    dest="structural_variants",
    default=False,
    help=
    '''Whether structural variants should be called and benchmarked with sniffles and truvari.
    ''')

    app.add_argument(
    "-trc", 
    "--tandem_repeat_calling",
    action="store_true",
    dest="tandem_repeats",
    default=False,
    help=
    '''Whether tandem repeats should be called and benchmarked with trgt and truvari.
    ''')

    app.add_argument(
    "--cmrg", 
    action="store_true", 
    dest="cmrg", 
    default=False,
    help=
    '''Additionally benchmark against the GIAB Challenging Medically Relevant
    Genes (CMRG v1.00) small-variant and SV benchmarks.
    ''')

    # app.add_argument(
    # "-hm", 
    # "--hifiasm",
    # action="store_false",
    # dest="hifiasm",
    # default = True,
    # help=
    # '''Turn off hifiasm (On by default)
    # ''')

    # app.add_argument(
    # "-cm", 
    # "--compleasm",
    # action="store_false",
    # dest="compleasm",
    # default = True,
    # help=
    # '''Turn off compleasm (On by default)
    # ''')

    # app.add_argument(
    # "-m", 
    # "--merqury",
    # action="store_false",
    # dest="merqury",
    # default = True,
    # help=
    # '''Turn off merqury (On by default)
    # ''')

    app.add_argument(
    "-ol", 
    "--outlier",
    action="store",
    dest="outlier",
    help=
    '''In combine mode, which of the input read sets should be treated as an outlier 
    (meaning it is much smaller than the rest and can not be downsampled to each of the desired fractions).
    ''')

    app.add_argument(
    "-t", 
    "--target_base_coverage",
    action="store",
    dest="coverage",
    type = int,
    help=
    '''For combine and downsample: To which target base coverage (number of nucleotides) will be downsampled.
    By default, will take the smallest read set present in the provided seqkit out-file as downsample target.
    ''')

    app.add_argument(
    "-nre", 
    "--no_restrict",
    action="store_false",
    dest="restrict",
    default = True,
    help=
    '''For combine and downsample: If the smallest given input is an outlier (< 1/3 the number of sequenced nts compared to 
    the biggest available read set), still use it to infer target base coverage. If not set,
    the next bigger read set is instead taken until one is found that is > 1/3 number of sequenced nucleotides of the biggest set.  
    ''')

    app.add_argument(
    "-dr",
    "--dry_run",
    action="store_true",
    default=False,
    help=
    """If set, a snakemake dry run with default parameters (--keep-going, --rerun-imcomplete, --the default DEF.yaml and prof/config.yaml files)
    will be started either after creation of the working directory, or, if the working directory already exists, within that working directory.
    """)

    app.add_argument(
    "-rs",
    "--run_snakemake",
    action="store_true",
    default=False,
    help=
    """If set, the snakemake run will not be executed as a dry-run, but run directly. Default is False.
    If one wants to execute the snakemake run, enable the -rs flag. If both --run_snakemake and --start_dry_run are set,
    only a dry run will be started.
    """)

    app.add_argument(
    "-sq",
    "--seqkit_file_path",
    action="store",
    dest="seqkit_path",
    help=
    """Path to seqkit file, based on which should be downsampled. (Will parse it and find the read set with the smallest output and that output)
    """)

    app.add_argument(
    "-sa",
    "--samples",
    action="store",
    dest="samples",
    help=
    """For combine: Take samples (must have same names as in seqkit file) to create combinations.
    Comma separated list, for example: readsA, readsB, readsC.
    """)

    app.add_argument(
    "-f",
    "--force_run",
    action="store_true",
    default=False,
    help=
    """If set, a polymeval command will be executed in the current wd. 
    Be careful, as any existing DEF file will be overwritten.
        """)
    
    app.add_argument(
    "--directory_name", 
    action="store",
    dest="directory_name",
    default = "polymeval_test",
    type=str, 
    help="""Name of the directory that polymeval pipeline is started in.
    Default is polymeval_test.
    """
    ) 

    app.add_argument(
    "--lo",
    "--local_run", 
    action="store_true",
    default=False,
    dest="local_run",
    help="""If no slurm scheduler is available, run the pipeline locally.
    """
    )  

    app.add_argument(
    "-rtt",
    "--rerun_triggers_mtime", 
    action="store_true",
    default=False,
    dest="rerun_trigger",
    help="""FOR DEVELOPMENT: If something in the polymeval code was changed, should reruns be done only on rules that have not yet produced proper output?
    (Snakemake --rerun-triggers mtime flag)
    """
    )  

    args = app.parse_args()
    return args


def main():
    args = argument_parser()

    ## Which mode to run
    if args.standard:
        logger.info("Polymeval will be run in standard mode (raw input reads)")
    elif args.downsample:
        logger.info("Polymeval will be run in downsample mode (reads will be downsampled first)")
    elif args.combine:
        logger.info("Polymeval will be run in combine mode (reads will be downsampled and then combined)")
    elif args.reference:
        logger.info("Polymeval will be run in reference mode (reads and assemblies will be compared to a reference assembly)")

    ## Create the YAML structure as a python dictionary.
    config = {
        "compleasm_db": "mollusca_odb12",
        "min_frac": 3,
        "kmc": False,
        "readstats": False,
        "hifieval": False,
        "hg_size": [],
        "remove_dups" : False,
        "colors": [],
        "reference_seq": "",
        "exons": [],
        "repeats": [],
        "combo_pairwise" : False,
        "gzipped": True,
        "structural_variants": False,
        "tandem_repeats": False,
        "wrapper_versions": {
            "meryl":"v9.4.2",
            "minimap":"v9.9.0",
            "bcftools":"v9.4.1",
            "tabix": "v9.4.1",
            "genomescope": "v9.4.2",
            "happy": "v7.0.0",
            "hifiasm": "v9.4.2",
            "samtools": "v9.4.2",
            "bbtools": "v9.16.0",
            "seqkit": "v9.4.2",
            "seqtk": "v7.0.0",
            "sniffles": "v9.18.0",
        }
    }

    ## Which assemblers are used in standard (and downsample mode)
    VALID = {"hifiasm", "flye", "lja", "verkko"}
    asms = [a.strip().lower() for a in args.assemblers.split(",") if a.strip()]
    bad = set(asms) - VALID
    if bad:
        logger.critical("Unknown assembler(s): %s. Choose from: %s", ", ".join(bad), ", ".join(sorted(VALID)))
        sys.exit(1)
    if "verkko" in asms:
        logger.warning("Verkko is designed for HiFi + ONT ultra-long reads. Running it on HiFi alone "
                    "is supported but costs far more compute than hifiasm for comparable "
                    "contiguity. See --verkko_no_correction.")
    config["assemblers"] = format_list(asms)
    config["verkko_extra"] = ""

    ## LJA binary wire-in
    if "lja" in asms:
        lja_bin = os.path.abspath(os.path.expanduser(args.lja_path)) if args.lja_path else "lja"
        if not (os.path.isfile(lja_bin) and os.access(lja_bin, os.X_OK)) and not shutil.which(lja_bin):
            logger.critical("LJA binary %r not found or not executable. Build it from source and pass "
                            "--lja_path, or drop 'lja' from --assembler.", lja_bin)
            sys.exit(1)
        config["lja_path"] = SingleQuotedScalarString(lja_bin)

    ## Fetching human benchmark data:
    config["cmrg"] = bool(args.cmrg)

    if args.fetch_benchmarks:
        if not args.benchmark_path:
            logger.critical("--fetch_benchmarks needs --benchmark_path.")
            sys.exit(1)
        os.makedirs(os.path.abspath(args.benchmark_path), exist_ok=True)
        config["benchmark_manifest"] = SingleQuotedScalarString(
            os.path.join(base_dir, "config", "benchmark.yaml"))
        config["benchmark_releases"] = format_list(
            [r.strip() for r in args.benchmark_releases.split(",") if r.strip()])

    ## Additional parameters:
    if (args.downsample or args.combine) and args.seqkit_path:
        seqkit_path = os.path.abspath(args.seqkit_path)
        config["path_to_seqkit"] =  SingleQuotedScalarString(seqkit_path)
    else:
        seqkit_path = None

    if (args.combine or args.downsample) and args.outlier:
        config["outlier"] = SingleQuotedScalarString(args.outlier)
    elif (args.combine or args.downsample) and not args.outlier:
        config["outlier"] = ""

    if (args.combine or args.downsample) and args.coverage:
        config["sample_base_target"] = int(args.coverage)
    elif (args.combine or args.downsample) and not args.coverage:
        config["sample_base_target"] = False

    if args.combine or args.downsample:
        config["restrict_downsampling"] = bool(args.restrict)

    ## Change default parameters
    compleasm_lib = resolve_compleasm_lib(args.compleasm_db_path)
    config["compleasm_db_path"] = SingleQuotedScalarString(compleasm_lib)
    if not args.compleasm_db_path:
        logger.info("No --compleasm_db_path given; using %s", compleasm_lib)
        logger.info("On a shared cluster, point --compleasm_db_path at a central ODB library, "
                    "or set POLYMEVAL_COMPLEASM_LIBS, to avoid one copy per user.")

    if args.compleasm_db:
        config["compleasm_db"] = SingleQuotedScalarString(args.compleasm_db)

    if args.kmc:
        config["kmc"] = True
        config["kmc_prep"] = get_kmc_prep_path()
    
    if args.hifieval:
        config["hifieval"] = True

    if args.readstats:
        config["readstats"] = True

    if args.remove_dups:
        config["remove_dups"] = True

    if args.reference_seq != "":
        reference_path = os.path.abspath(args.reference_seq)
        config["reference_seq"] = SingleQuotedScalarString(reference_path)
    else:
        reference_path = None

    if args.colors != "":
        color_path = os.path.abspath(args.colors)
        config["colors"] = SingleQuotedScalarString(color_path)
    else:
        color_path = None

    if args.hg_size:
        config["hg_size"] = args.hg_size

    if args.pairwise:
        config["combo_pairwise"] = True

    if args.rasusa_seed:
        config["seed"] = int(args.rasusa_seed)

    if args.benchmark_path:
        benchmark_path = os.path.abspath(args.benchmark_path)
        config["benchmark_path"] = SingleQuotedScalarString(benchmark_path)
    else:
        benchmark_path = None

    if args.structural_variants:
        config["structural_variants"] = True

    if args.tandem_repeats:
        config["tandem_repeats"] = True

    ## Set up directory;
    READS_SUBDIR = "raw_reads"
    SUFFIXES = tuple(sorted((".dup.fastq.gz", ".dup.fastq", ".dup.fq.gz", ".dup.fq") if config["remove_dups"]
                        else (".fastq.gz", ".fastq", ".fq.gz", ".fq"), key=len, reverse=True))
    GZ  = (".fastq.gz", ".fq.gz")

    def link_and_discover(src, work_dir, suffixes, wanted=None, reference_run=False, skip_bases_ending=None):
        dest = os.path.join(work_dir, READS_SUBDIR)
        found = link_reads(src, dest, suffixes, down_list=wanted,
                        skip_bases_ending=skip_bases_ending)

        if not found:
            logger.critical("No files ending in %s found in %s.", " or ".join(suffixes), src)
            sys.exit(1)

        if wanted is not None:
            missing = wanted - {b for b, _ in found}
            if missing:
                logger.critical("Requested sample(s) not found in %s: %s",
                                src, ", ".join(sorted(missing)))
                sys.exit(1)

        return dest, found

    def resolve_gzipped(found, dest):
        """All inputs must be compressed, or none. Returns True/False."""
        compressed = {suf.endswith(".gz") for _, suf in found}
        if len(compressed) > 1:
            logger.critical("Directory %s mixes compressed and uncompressed reads. polymeval needs "
                            "all inputs in the same form — gzip the plain files and rerun.", dest)
            sys.exit(1)
        return compressed.pop()

    def basenames(found):
        return format_list([b for b, _ in found])

    work_dir = os.path.join(os.getcwd(), args.directory_name)
    os.makedirs(work_dir, exist_ok=True)
    wanted = {x.strip() for x in args.samples.split(",") if x.strip()} if args.samples else None

    if args.fetch_benchmarks:
        pass

    ## GZ as suffix here because rules only allow gz files and no deduplication
    elif args.combine or args.downsample:
        readset_dict, downsample_samples, removed_samples, downsample_nucs = \
            get_downsample_rates.read_seq_stats(args.seqkit_path,
                                                restrict=config["restrict_downsampling"],
                                                min_frac=config["min_frac"])
        if args.combine:
            if not wanted:
                logger.critical("--combine requires --samples."); sys.exit(1)
            path_for_link_rds, found = link_and_discover(args.in_reads, work_dir, GZ, wanted=wanted)
            config["samples"] = format_list(sorted(wanted))
        else:
            path_for_link_rds, found = link_and_discover(args.in_reads, work_dir, GZ,
                                                        wanted=wanted or set(downsample_samples))
            config["samples"] = basenames(found)

    elif args.reference:
        path_for_link_rds, found = link_and_discover(args.in_reads, work_dir, GZ,
                                                    skip_bases_ending=".dup")
        path_for_link_asm = os.path.join(work_dir, "assemblies")
        asm_found = link_assemblies(args.in_assemblies, path_for_link_asm)

        reads = {b for b, _ in found}
        asms  = {b for b, _ in asm_found}
        if reads != asms:
            logger.critical("Reads and assemblies do not match.\n  only in reads: %s\n"
                            "  only in assemblies: %s",
                            ", ".join(sorted(reads - asms)) or "-",
                            ", ".join(sorted(asms - reads)) or "-")
            sys.exit(1)
        config["samples"] = basenames(found)

    # standard and variant-calling
    else:                                      
        path_for_link_rds, found = link_and_discover(args.in_reads, work_dir, SUFFIXES, wanted=wanted, skip_bases_ending=None if config["remove_dups"] else ".dup")
        config["gzipped"] = resolve_gzipped(found, path_for_link_rds)
        config["samples"] = basenames(found)

    ## Which snakefile to use:
    if args.standard:
        snakefile = "Snakefile_standard"

    elif args.downsample:
        ## Change so that user can combine downsampling target
        if not args.seqkit_path:
            logger.critical("If downsampling mode should be run, please specify a seqkit in-file")
            sys.exit(1)
        snakefile = "Snakefile_downsample"

    elif args.combine:
        if not args.samples:
            logger.critical("If combine mode should be run, please specify a which samples should be used for combinations")
            sys.exit(1)
        snakefile = "Snakefile_combine"
    
    elif args.reference:
        if not args.in_reads or not args.in_assemblies:
            logger.critical("No reads or assemblies missing! Please specify both with --input_reads and --input_assemblies")
            sys.exit(1)
        if args.reference_seq == "":
            logger.critical("A reference genome sequence has to be provided to --reference_seq in order to run the analysis.")
            sys.exit(1)
        snakefile = "Snakefile_reference"
    
    elif args.variant_calling:
        if not args.benchmark_path:
            logger.critical("No path to benchmark files given. You can set the path to benchmark files with --benchmark_path")
            sys.exit(1)
        elif not os.path.exists(args.benchmark_path):
            logger.critical("Path to benchmark files given by --benchmark_path does not exist. Check that the path was spelled correctly and exists")
            sys.exit(1)
        snakefile = "Snakefile_human_vcf"

    elif args.fetch_benchmarks:
        snakefile = "Snakefile_fetch"


    ## Add DEF file to created directory
    # Write to the specified YAML file
    yaml_name = 'DEF.yaml'
    
    ## Get working directory for the current run to put DEF file in:
    path_for_DEF = os.path.join(os.getcwd(), args.directory_name, yaml_name)    
    
    ## Check whether a DEF.yaml file already exists. If it does, first sys.exit, unless the --force-run is enabled
    ## Whether DEF file has been updated
    update_DEF = False

    if os.path.exists(path_for_DEF):
        current_yaml = YAML()
        with open(path_for_DEF, 'r') as file:
            current_def_file = current_yaml.load(file)
        if current_def_file == config:
            logger.info("This command has been run before. The DEF.yaml file will be overwritten, but its contents will not change.")

        elif current_def_file != config and args.force_run == False:
            logger.critical("DEF.yaml already exists in current working directory :{}, its content differs from your input commands. If you want to force a run in this working directory, enable -f/--force_run".format(path_for_DEF))
            sys.exit(1)

        elif current_def_file != config and args.force_run:
            logger.warning("""
    WARNING: There is an existing polymeval directory in the specified working directory: {}. 
    The commands with which it was initialized differ from the currently invoked commands, 
    but the run will be FORCED with the currently enabled commands.
    """.format(path_for_DEF.split('DEF.yaml')[0]))
            update_DEF = True

    with open(path_for_DEF, "w") as yaml_file:
        yaml = YAML()
        yaml.boolean_representation = ['False', 'True']
        yaml.default_flow_style = False
        #yaml.indent(mapping=2, sequence=4, offset=2)
        yaml.preserve_quotes = True
        yaml.dump(config, yaml_file)

    logger.info(f"YAML configuration written to DEF.yaml")

    ## Whether apptainer should be used:
    if args.variant_calling or args.reference:
        use_apptainer = True
    else:
        use_apptainer = False

    ## get apptainer arguments
    if args.variant_calling:
        apptainer_args = f"-B ./alignments:/input -B {benchmark_path}:/reference -B ./variants:/output"
    elif args.reference:
        apptainer_args = get_apptainer_bind_args([
            args.reference_seq,
        ])
    else:
        apptainer_args = None

    if args.dry_run:
        run_snakemake(snake_file = snakefile, 
                        config_file = None, 
                        directory_name = args.directory_name, 
                        dryrun = args.dry_run, 
                        snake_default = True,
                        conda_path = "",
                        rerun_triggers = args.rerun_trigger, 
                        updated_rule = update_DEF,
                        local_run = args.local_run,
                        use_apptainer = use_apptainer,
                        apptainer_args = apptainer_args)

    elif args.run_snakemake:
        run_snakemake(snake_file = snakefile, 
                        config_file = None, 
                        directory_name = args.directory_name, 
                        dryrun = not args.run_snakemake, 
                        snake_default = True,
                        conda_path = "",
                        rerun_triggers = args.rerun_trigger,
                        updated_rule = update_DEF,
                        local_run = args.local_run,
                        use_apptainer = use_apptainer,
                        apptainer_args = apptainer_args)


if __name__ == "__main__":
    main()