import argparse
import sys
import os
from ruamel.yaml import load_all, YAML, comments
from ruamel.yaml.scalarstring import SingleQuotedScalarString, DoubleQuotedScalarString
import subprocess
import get_downsample_rates

## Absolute paths on file system
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "."))

def symlink_all_rds(src_path, dst_path, down_list = []): 
    os.makedirs(dst_path, exist_ok=True)   

    for rds in os.listdir(src_path):
        if down_list == []:
            down_list = os.listdir(src_path)

        if rds.endswith("fastq.gz") and (rds.split('.')[0] in down_list or rds in down_list):
            src_file = os.path.join(os.path.realpath(src_path), rds)
            dst_file = os.path.join(dst_path, rds)
            if not os.path.islink(dst_file):
                os.symlink(src_file, dst_file)

def symlink_all_asm(src_path, dst_path): 
    os.makedirs(dst_path, exist_ok=True)   

    for asm in os.listdir(src_path):
        if asm.endswith("fa"):
            src_file = os.path.join(os.path.realpath(src_path), asm)
            dst_file = os.path.join(dst_path, asm)
            if not os.path.islink(dst_file):
                os.symlink(src_file, dst_file)

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
                  local_run = False):
    
    ## move to the working directory
    current_wd =  os.getcwd()
    if working_dir == current_wd:
        snakemake_dir = os.path.join(current_wd, directory_name)
    else:
        snakemake_dir = os.path.join(working_dir, directory_name)

    if not os.path.exists(snakemake_dir):
        os.makedirs(snakemake_dir)
    
    #snake_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "."))
    #config_file = 

    os.chdir(snakemake_dir)
    
    ## Run snakemake from the command line
    cmd = ['snakemake']

    snake_file = get_snakefile_path(snake_file)
    cmd += ['--snakefile', snake_file]
    
    if not local_run:
        config_file = get_cluster_configfile_path()
        cmd += ['--profile', config_file]
    else:
        ## Add these arguments to CLI so users can set them
        cmd += ['--cores', '4', '--max-threads', '4', '--resources', 'mem_mb=30000']

    cmd += ['--directory', snakemake_dir]

    ## no lock to allow reruns
    cmd += ['--nolock'] 

    if use_conda and not conda_path:
        conda_path = os.path.join(base_dir, ".snakemake/conda")
        cmd += ['--use-conda', '--conda-prefix', conda_path]
    elif use_conda and conda_path:
        cmd += ['--use-conda', '--conda-prefix', conda_path]
    elif not use_conda and conda_path:
        sys.exit("--use-conda was not specified, but a reference conda path given. Since use_conda is activated by default, do not activate it if you also pass a conda path")
    elif not use_conda and not conda_path:
        print("WARNING: --use-conda is deactivated and no conda path is found. Most likely the pipeline will fail, unless all tools are installed on the user's machine and availabe in $PATH")

    if dryrun:
        cmd.append('--dry-run')
    
    if rerun_triggers:
        cmd += ['--rerun-triggers', 'mtime']    
    
    if snake_default:
        default_snakemake_args = ["--rerun-incomplete", "--keep-going"]
        cmd += default_snakemake_args

    print_cmd = " ".join(cmd)
    print("The following command will be run for the polymeval pipeline: {}".format(print_cmd))

    ## Run pipeline proper
    try:
        result1 = subprocess.run(cmd, check=True)
        print(result1.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")

DESCRIPTION = '''
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
    "-ia", 
    "--input_assemblies",
    action="store",
    dest="in_assemblies",
    type=str,
    help=
    '''For 'reference' mode: Path to directory with input assemblies. Naming as with input reads, and names have to match.
    Normally, these assemblies will have been created with --standard, --combine or --downsample and the naming therefore
    will work out naturally. For the best experience, pandepth should be installed on the system and accessible in path.
    The path to pandepth can also be passed through --pandepth_path.
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
    Must be present in the specified path, otherwise will be downloaded automatically 
    (It is highly advised to have one accessible odb library stored centrally on the HPC to
    circumvent repeated costly downloads), see --compleasm_db_path
    ''')

    app.add_argument(
    "-dbp", 
    "--compleasm_db_path",
    action="store",
    dest="compleasm_db_path",
    default="/sgn/software/orthodb/odb12_latest/",
    type=str,
    help=
    '''Path to accessible odb library stored centrally on the HPC. If not provided, 
    compleasm will download libraries on the fly, which takes a lot of time and disk space.
    ''')

    app.add_argument(
    "-co", 
    "--colors",
    action="store",
    dest="colors",
    default="",
    type=str,
    help=
    '''Absolute path to a tsv file listing sample names in column 1 and 
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

    app.add_argument(
    "-pa", 
    "--pandepth",
    action="store_true",
    dest="pandepth",
    help=
    '''Turn on Pandepth depth mapping. Most reference mode analyses will need this.
    ''')

    app.add_argument(
    "-pp", 
    "--pandepth_path",
    action="store",
    dest="pandepth_path",
    default = "",
    help=
    '''In case pandepth is installed but not in the user's $PATH, provide absolute path to pandepth.
    ''')

    app.add_argument(
    "-km", 
    "--kmer_length",
    action="store",
    dest="kmer_length",
    default=25,
    help=
    '''K-mer length used by meryl and KMC (if enabled). Default is 25.
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
    "-re", 
    "--restrict",
    action="store_true",
    dest="restrict",
    default = True,
    help=
    '''For combine and downsample: If the smallest given input is an outlier (< 1/3 the number of sequenced nts compared to 
    the biggest available read set), whether it should still be used to infer target base coverage. If -r is set, 
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

    args = app.parse_args()
    return args


def main():
    args = argument_parser()

    ## Which mode to run
    if args.standard:
        print("Polymeval will be run in standard mode (raw input reads)")
    elif args.downsample:
        print("Polymeval will be run in downsample mode (reads will be downsampled first)")
    elif args.combine:
        print("Polymeval will be run in combine mode (reads will be downsampled and then combined)")
    elif args.combine:
        print("Polymeval will be run in reference mode (reads and assemblies will be compared to a reference assembly)")

    ## Create the YAML structure as a python dictionary.
    config = {
        "compleasm_db": "mollusca_odb12",
        "min_frac": 3,
        "kmc": False,
        "readstats": False,
        "hifieval": False,
        "remove_dups" : False,
        "colors": [],
        "pandepth": False,
        "reference_seq": "",
        "exons": [],
        "repeats": [],
        "pandepth_path": [],
    }

    ## Additional parameters:
    if (args.downsample or args.combine) and args.seqkit_path:
        config["path_to_seqkit"] =  os.path.abspath(args.seqkit_path)

    if (args.combine or args.downsample) and args.outlier:
        config["outlier"] = SingleQuotedScalarString(args.outlier)
    elif (args.combine or args.downsample) and not args.outlier:
        config["outlier"] = ""

    if (args.combine or args.downsample) and args.coverage:
        config["sample_base_target"] = int(args.coverage)
    elif (args.combine or args.downsample) and not args.coverage:
        config["sample_base_target"] = ""

    if (args.combine or args.downsample) and args.restrict:
        config["restrict_downsampling"] = True
    elif (args.combine or args.downsample) and not args.coverage:
        config["sample_base_target"] = False

    ## Change default parameters
    if args.compleasm_db:
        config["compleasm_db"] = SingleQuotedScalarString(args.compleasm_db)

    if args.compleasm_db_path:
        config["compleasm_db_path"] = SingleQuotedScalarString(args.compleasm_db_path)

    if args.kmc:
        config["kmc"] = True
        config["kmc_prep"] = get_kmc_prep_path()
    
    if args.hifieval:
        config["hifieval"] = True

    if args.readstats:
        config["readstats"] = True

    if args.pandepth:
        config["pandepth"] = True

    if args.pandepth_path != "":
        config["pandepth_path"] = SingleQuotedScalarString(args.pandepth_path)

    if args.reference_seq != "":
        config["reference_seq"] = SingleQuotedScalarString(args.reference_seq)

    if args.colors != "":
        config["colors"] = SingleQuotedScalarString(args.colors)

    ## Set up directory;
    dest_path = os.path.join(os.getcwd(), args.directory_name)
    if not os.path.exists(dest_path):
        os.makedirs(dest_path, exist_ok=True)  

    if not args.standard and not args.reference:
        readset_dict, downsample_samples, removed_samples, downsample_nucs = get_downsample_rates.read_seq_stats(args.seqkit_path, restrict = config["restrict_downsampling"], min_frac = config["min_frac"])

        if (args.combine and args.samples):
            path_for_link_rds = os.path.join(os.getcwd(), args.directory_name, "raw_reads")
            samples = format_list(args.samples.split(','))
            symlink_all_rds(args.in_reads, path_for_link_rds, samples)
            in_reads = os.listdir(path_for_link_rds)
            in_reads_list = [f.replace('.fastq.gz','') for f in in_reads if (os.path.islink(os.path.join(path_for_link_rds, f)) or os.path.isfile(os.path.join(path_for_link_rds, f))) and f.endswith(".fastq.gz") and f.replace(".fastq.gz", "") in samples]
            config["samples"] =  samples

        elif (args.downsample and not args.samples):
            path_for_link_rds = os.path.join(os.getcwd(), args.directory_name, "raw_reads")
            symlink_all_rds(args.in_reads, path_for_link_rds, downsample_samples)
            in_reads = os.listdir(path_for_link_rds)
            in_reads_list = [f.replace('.fastq.gz','') for f in in_reads if (os.path.islink(os.path.join(path_for_link_rds, f)) or os.path.isfile(os.path.join(path_for_link_rds, f))) and f.endswith(".fastq.gz") and f.replace(".fastq.gz", "") in downsample_samples]
            config["samples"] =  format_list(in_reads_list)
        
    elif args.reference:
        path_for_link_rds = os.path.join(os.getcwd(), args.directory_name, "raw_reads")
        symlink_all_rds(args.in_reads, path_for_link_rds, [])
        path_for_link_asm = os.path.join(os.getcwd(), args.directory_name, "assemblies")
        symlink_all_asm(args.in_assemblies, path_for_link_asm)
        in_reads = os.listdir(path_for_link_rds)
        in_reads_list = [f.replace('.fastq.gz','') for f in in_reads if (os.path.islink(os.path.join(path_for_link_rds, f)) or os.path.isfile(os.path.join(path_for_link_rds, f))) and f.endswith(".fastq.gz")]
        in_assemblies = os.listdir(path_for_link_asm)
        in_assemblies_list = [f.replace('.fa','') for f in in_assemblies if (os.path.islink(os.path.join(path_for_link_asm, f)) or os.path.isfile(os.path.join(path_for_link_asm, f))) and f.endswith(".fa")]
        if sorted(in_reads_list) == sorted(in_assemblies_list):
            config["samples"] =  format_list(in_reads_list)
        else:
            sys.exit("Number or names of assemblies and raw reads differ. Please check that they have the same base names and that there is the same number of them present in the respective directories")
    
    else:
        if args.samples:
            path_for_link_rds = os.path.join(os.getcwd(), args.directory_name, "raw_reads")
            samples = format_list(args.samples.split(','))
            symlink_all_rds(args.in_reads, path_for_link_rds, samples)
            in_reads = os.listdir(path_for_link_rds)
            in_reads_list = [f.replace('.fastq.gz','') for f in in_reads if (os.path.islink(os.path.join(path_for_link_rds, f)) or os.path.isfile(os.path.join(path_for_link_rds, f))) and f.endswith(".fastq.gz") and f.replace(".fastq.gz", "") in samples]
            config["samples"] =  format_list(in_reads_list)
        else:
            path_for_link_rds = os.path.join(os.getcwd(), args.directory_name, "raw_reads")
            symlink_all_rds(args.in_reads, path_for_link_rds, [])
            in_reads = os.listdir(path_for_link_rds)
            in_reads_list = [f.replace('.fastq.gz','') for f in in_reads if (os.path.islink(os.path.join(path_for_link_rds, f)) or os.path.isfile(os.path.join(path_for_link_rds, f))) and f.endswith(".fastq.gz")]
            config["samples"] =  format_list(in_reads_list)


    ## Which snakefile to use:
    if args.standard:
        snakefile = "Snakefile_standard"

    elif args.downsample:
        ## Change so that user can combine downsampling target
        if not args.seqkit_path:
            sys.exit("If downsampling mode should be run, please specify a seqkit in-file")
        snakefile = "Snakefile_downsample"

    elif args.combine:
        if not args.samples:
            sys.exit("If combine mode should be run, please specify a which samples should be used for combinations")
        snakefile = "Snakefile_combine"
    
    elif args.reference:
        if not args.in_reads or not args.in_assemblies:
            sys.exit("No reads or assemblies missing! Please specify both with --input_reads and --input_assemblies")
        if args.reference_seq == "":
            sys.exit("A reference genome sequence has to be provided to --reference_seq in order to run the analysis.")
        if not args.pandepth:
            print("Warning: Pandepth mode is not turned on, most of the analysis cannot be run. Install pandepth for the best experience")
        if args.pandepth_path == "":
            print("Warning: Pandepth path not given. If pandepth is not installed or accessible, all operations involving pandepth will fail.")
        snakefile = "Snakefile_reference"


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
        #print(current_def_file == config)
            print("Info: This command has been run before. The DEF.yaml file will be overwritten, but its contents will not change.")

        elif current_def_file != config and args.force_run == False:
        #print(current_def_file == config)
            sys.exit("DEF.yaml already exists in current working directory :{}, its content differs from your input commands. If you want to force a run in this working directory, enable -f/--force_run".format(path_for_DEF))
    
        elif current_def_file != config and args.force_run:
            print("""
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

    print(f"YAML configuration written to DEF.yaml")

    if args.dry_run:
        run_snakemake(snake_file = snakefile, 
                        config_file = None, 
                        directory_name = args.directory_name, 
                        dryrun = args.dry_run, 
                        snake_default = True,
                        conda_path = "",
                        #rerun_triggers = args.rerun_trigger, 
                        updated_rule = update_DEF,
                        local_run = args.local_run)

    elif args.run_snakemake:
        run_snakemake(snake_file = snakefile, 
                        config_file = None, 
                        directory_name = args.directory_name, 
                        dryrun = not args.run_snakemake, 
                        snake_default = True,
                        conda_path = "",
                        #rerun_triggers = args.rerun_trigger,
                        updated_rule = update_DEF,
                        local_run = args.local_run)


if __name__ == "__main__":
    main()