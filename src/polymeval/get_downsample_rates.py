import argparse
import logging
from collections import defaultdict
from itertools import combinations
import sys

## Logging
logger = logging.getLogger("downsample_rates")
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("[%(levelname)-8s] %(message)s"))
logger.addHandler(handler)
logger.propagate = False  # ← critical: stops messages bubbling up to root logger

DESCRIPTION = ""

def parse_args():
    """Parse CMD args."""

    app = argparse.ArgumentParser(description=DESCRIPTION)
    
    app.add_argument(
    "-s", 
    "--seqkit_stats_file",
    action="store",
    dest="seqkit",
    help="Path to seqkit stats file")

    app.add_argument(
    "-c", 
    "--combinations",
    action="store",
    dest="combinations",
    help="List of read sets that should be combined")

    app.add_argument(
    "-r", 
    "--restrict",
    action="store_true",
    dest="restrict",
    help="Whether to restrict downsampling to a minimum, for example, if the minimum number of input nts is less than 1/2 the maximum, take the next higher number as target")

    app.add_argument(
    "--pairwise-single",
    action="store_true",
    dest="pairwise_single",
    help="Enable mode to use up to 10 inputs, downsample to 'single' level, but form pairwise combinations")

    app.add_argument(
    "-m", 
    "--min_frac",
    action="store",
    dest="min_frac",
    default=5,
    help="Minimum fraction of max sample for a sample to be included")

    app.add_argument(
    "-bt", 
    "--base_target",
    action="store",
    dest="base_target",
    help="sample base target (to what amount should be downsampled)")

    app.add_argument(
    "-ol", 
    "--outlier",
    action="store",
    dest="outlier",
    help="Sample(s) that might now match the sample base target")
       
    args = app.parse_args()

    return args

## Global variables
## names for input files in combination trials and number of files to be combined for those names
suffix_dict = {1 : "single", 2 : "half", 3 : "third", 4 : "fourth", 5: "fifth", 6: "sixth"}
combination_group_sizes = {"single" : 1, "half": 2, "third": 3, "fourth": 4, "fifth": 5, "sixth": 6}

def read_seq_stats(path, restrict, min_frac, samples=None):
    readset_dict = {}
    all_sample_list = []
    with open(path) as s:
        for line in s:
            if line.startswith("file"):
                continue
            else:
                stat_line = line.strip().split('\t')
                if "/" in stat_line[0]:
                    name = stat_line[0].split('/')[1].split('.')[0]
                else:
                    name = stat_line[0].split('.')[0]
                nucs = int(stat_line[4])
                readset_dict[name] = nucs
                all_sample_list.append(name)

    if samples is not None:
        missing = [s for s in samples if s not in readset_dict]
        if missing:
            logger.critical("Requested sample(s) not in %s: %s\n  present: %s",
                            path, ", ".join(sorted(missing)), ", ".join(sorted(readset_dict)))
            sys.exit(1)

        ignored = [s for s in readset_dict if s not in samples]
        if ignored:
            logger.info("Using the %d requested read set(s); ignoring %d other row(s) in %s: %s",
                        len(samples), len(ignored), path, ", ".join(sorted(ignored)))

        readset_dict = {s: readset_dict[s] for s in samples}

        ## An explicit list is a choice, not a candidate pool. Report the spread, drop nothing.
        lo_k = min(readset_dict, key=readset_dict.get)
        hi_k = max(readset_dict, key=readset_dict.get)
        if restrict and readset_dict[lo_k] < readset_dict[hi_k] / min_frac:
            logger.warning("%s (%s nt) is below 1/%s of %s (%s nt). It was named explicitly, so it "
                           "is kept and sets the target; every other read set is reduced to it. "
                           "Drop it from --samples, or pass --target_bases, if that is not "
                           "what you want.",
                           lo_k, f"{readset_dict[lo_k]:,}", min_frac,
                           hi_k, f"{readset_dict[hi_k]:,}")
        restrict = False


    ## Guard against sinlge-read entries:
    if len(readset_dict) < 2:
        logger.critical("Downsampling needs at least two read sets; %s describes %d.",
                        path, len(readset_dict))
        sys.exit(1)

    min_read_set = adjust_min_max(readset_dict, restrict, min_frac)
    min_read_set_size = readset_dict[min_read_set]

    logger.info("Downsample target: %s nt (smallest of %d read set(s): %s)",
            f"{min_read_set_size:,}", len(readset_dict), min_read_set)

    samples_out = [r for r in readset_dict if r != min_read_set]
    removed = [r for r in all_sample_list if r not in readset_dict]
    return readset_dict, samples_out, removed, min_read_set_size
    

def adjust_min_max(readset_dict, restrict, min_frac):
    while restrict and len(readset_dict) > 2:
        min_rs = min(readset_dict, key=readset_dict.get)
        max_rs = max(readset_dict, key=readset_dict.get)
        if readset_dict[min_rs] >= readset_dict[max_rs] / min_frac:
            break
        logger.warning("Read set %s (%d nt) is below 1/%s of the largest, %s (%d nt). "
                       "Excluding it; the next smallest becomes the downsampling target.",
                       min_rs, readset_dict[min_rs], min_frac, max_rs, readset_dict[max_rs])
        readset_dict.pop(min_rs)

    min_read_set = min(readset_dict, key=readset_dict.get)
    max_read_set = max(readset_dict, key=readset_dict.get)

    if restrict and readset_dict[min_read_set] < readset_dict[max_read_set] / min_frac:
        logger.warning("Only %d read sets remain, so %s is kept as the downsampling target "
                       "despite being below 1/%s of %s. Every other sample will be reduced "
                       "to %d nt.", len(readset_dict), min_read_set, min_frac, max_read_set,
                       readset_dict[min_read_set])
    return min_read_set

def create_combination_downsamples(readset_dict, read_sets, outlier, minimum, pairwise_single=False):
    number_combos = len(read_sets)

    if pairwise_single:
        logger.info("Only pairwise combinations will be run")

    if outlier:
        logger.info("{} is defined as an outlier read set".format(outlier))

    if minimum:
        logger.info("minimum downsampling target is {}".format(minimum))

    ## Cases than cannot be handled
    max_combos = 15 if pairwise_single else 6
    if number_combos > max_combos:
        logger.critical(f"You can at most combine {max_combos} read sets")
        sys.exit(1)

    if number_combos < 2:
        logger.critical("There are fewer than one read sets available, please inset at least two read sets to be combined/downsampled")
        sys.exit(1)

    ## Tell the user about potentially exploding combinations
    warn_if_large(number_combos, pairwise_single)

    new_readset_dict = {key: readset_dict[key] for key in read_sets}
    new_min_read_set = min(new_readset_dict, key=new_readset_dict.get)

    if minimum:
        target_total  = minimum
        target_source = "--target_bases"
        #new_min_read_set_size = minimum

    elif outlier and new_min_read_set == outlier:
        tmp_readset_dict = {key: readset_dict[key] for key in read_sets if key != outlier}
        new_min_read_set = min(tmp_readset_dict, key=tmp_readset_dict.get)
        target_total   = tmp_readset_dict[new_min_read_set]
        target_source = f"smallest non-outlier: {new_min_read_set}"
        
    else:
        target_total   = new_readset_dict[new_min_read_set]
        target_source = f"smallest selected read set: {new_min_read_set}"
    
    ## Pairwise standard mode is that half the target cov is used
    files_per_combo = 2 if pairwise_single else 1
    base_per_file   = target_total // files_per_combo

    downsample_nucs = []
    downsample_dict = {}
    downsample_read_dict = defaultdict(list)
    downsample_read_frac = defaultdict(list)

    ## Dynamic generation of levels
    nuc_list = []
    turned_nuc_dict = {}

    ## Standard mode: divisors go from 1 to n
    ## Pairwise mode: Divisors go from 1 to 2 (1/2)

    divisors = range(1, number_combos + 1)
    if pairwise_single:
        divisors = range(1,2)
    
    for div in divisors:
        if div == 1:
            n_nucs = base_per_file
        else:
            n_nucs = int(base_per_file / div)
        
        suffix = suffix_dict[div]
        nuc_list.append(n_nucs)
        downsample_dict[suffix] = n_nucs
        turned_nuc_dict[n_nucs] = suffix
        downsample_nucs.append(n_nucs)

    ## downsample operations:
    skipped = defaultdict(list)
    for read_set, size in new_readset_dict.items():
        for number in nuc_list:
            if size < number:
                skipped[read_set].append(turned_nuc_dict[number])
                continue
            downsample_read_dict[read_set].append(number)
            downsample_read_frac[read_set].append(turned_nuc_dict[number])

    ## Plan for logging:
    plan = dict(
        readset_dict    = dict(new_readset_dict),
        target_total    = target_total,
        target_source   = target_source,
        downsample_dict = dict(downsample_dict),
        pairwise_single = pairwise_single,
        skipped         = dict(skipped),
    )

    return downsample_read_dict, downsample_read_frac, downsample_nucs, downsample_dict, plan


def form_combinations(downsample_dict, downsample_read_frac, pairwise_single=False):
    file_names = [key for key in set(downsample_dict.keys())]
    
    combination_ranges = list(range(1, len(file_names)+1))

    ## In pairwise combos, only use single-suffix
    if pairwise_single:
        suffixes = ["single"]
    else:
        suffixes = [suffix_dict[number] for number in combination_ranges]

    all_combos = {}

    for sf in suffixes:
        #subset_files = [f"{f}.downsampled.{sf}.fastq.gz" for f in file_names]
        subset_files = [f"{f}-{sf}" for f in file_names if sf in downsample_read_frac[f]]
        n = combination_group_sizes[sf]

        ## For pairwise combinations
        if pairwise_single and sf == "single":
            n = 2

        if n == 1:
            for subset_file in sorted(subset_files):
                combo =  subset_file + ".fastq.gz"
                outname = "PLUS".join(subset_file.split("-"))
                all_combos[combo] = outname
        else:
            for combo in combinations(sorted(subset_files), n):
                combo_fastq = " ".join([fastq + '.fastq.gz' for fastq in combo])
                outname = "PLUS".join(combo)
                all_combos[combo_fastq] = outname

    # print(f"Generated {len(all_combos)} combinations:")
    # for c in all_combos:
    #     print(c)
    
    return(all_combos)

## Warning message for large combos:
def warn_if_large(number_combos, pairwise_single):
    """Say what the matrix will cost before building it."""
    if pairwise_single:
        n_combos = number_combos * (number_combos - 1) // 2
        rasusa   = number_combos
        kind     = "pairwise"
    else:
        n_combos = 2 ** number_combos - 1
        rasusa   = number_combos ** 2
        kind     = "all-vs-all"

    if n_combos > 50:
        logger.warning("%d read sets give %d %s combinations: %d assemblies, %d compleasm "
                       "runs, ~%d jobs. Check this is intended.",
                       number_combos, n_combos, kind, n_combos, n_combos,
                       7 * n_combos + rasusa + 3)

## Comprehensive Logging 
def log_plan(plan, n_combos):
    readset_dict    = plan["readset_dict"]
    downsample_dict = plan["downsample_dict"]
    pairwise_single = plan["pairwise_single"]
    target_total    = plan["target_total"]
    target_source   = plan["target_source"]
    skipped         = plan["skipped"]
    files_per = {s: (2 if pairwise_single else combination_group_sizes[s])
                 for s in downsample_dict}

    logger.info("Downsampling plan (%s)", "pairwise" if pairwise_single else "combine")
    logger.info("  read sets (%d):", len(readset_dict))
    for k, v in sorted(readset_dict.items(), key=lambda kv: -kv[1]):
        logger.info("    %-28s %15s nt", k, f"{v:,}")
    logger.info("  target per combination: %s nt (%s)", f"{target_total:,}", target_source)
    logger.info("  levels:")
    for suffix in sorted(downsample_dict, key=lambda s: combination_group_sizes[s]):
        per, n = downsample_dict[suffix], files_per[suffix]
        logger.info("    %-7s %d file(s) x %15s nt = %15s nt",
                    suffix, n, f"{per:,}", f"{per * n:,}")
    logger.info("  %d combination(s) will be built", n_combos)

    for read_set, levels in sorted(skipped.items()):
        logger.warning("  %s (%s nt) is too small for level(s) %s and will not appear in them.",
                       read_set, f"{readset_dict[read_set]:,}", ", ".join(levels))

## Not used for now
def main():
    args = parse_args()
    seqkit = args.seqkit
    restrict = args.restrict
    if args.combinations:
        combinations = args.combinations.split(',')
    else:
        logger.info("No combinations are to be done, instead just downsample all reads to the mininum input")

    readset_dict, samples, removed_samples, downsample_nucs = read_seq_stats(seqkit, restrict)
    
    if args.combinations:
        downsample_dict, new_min_read_set, new_min_read_set_size, plan = create_combination_downsamples(readset_dict, combinations)
        form_combinations(downsample_dict)

    log_plan

if __name__ == '__main__':
    main()




