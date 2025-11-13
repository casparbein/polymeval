import argparse
from collections import defaultdict
from itertools import combinations
import sys

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
       
    args = app.parse_args()

    return args

## Global variables
## names for input files in combination trials and number of files to be combined for those names
suffix_dict = {1 : "single", 2 : "half", 3 : "third", 4 : "fourth", 5: "fifth"}
combination_group_sizes = {"single" : 1, "half": 2, "third": 3, "fourth": 4, "fifth": 5}

def adjust_min_max(readset_dict, restrict, min_frac):
    min_read_set = min(readset_dict, key=readset_dict.get)
    max_read_set = max(readset_dict, key=readset_dict.get)

    if restrict and readset_dict[min_read_set] < readset_dict[max_read_set] / min_frac:
        print("WARNING: Minimum read set {} is smaller ({} nts) than 1/{} of the maximum {} ({} nts), setting minimum to the next suitable minimum".format(min_read_set, readset_dict[min_read_set], min_frac, max_read_set, readset_dict[max_read_set]))

        readset_dict.pop(min_read_set)

        return adjust_min_max(readset_dict, restrict,min_frac)
    
    else:
        return min_read_set


def read_seq_stats(path, restrict,min_frac):
    readset_dict = {}
    all_sample_list = []
    with open(path) as s:
        for line in s:
            if line.startswith("file"):
                continue
            else:
                stat_line = line.strip().split('\t')
                #print("/" in stat_line[0])
                if "/" in stat_line[0]:
                    name = stat_line[0].split('/')[1].split('.')[0]
                else:
                    name = stat_line[0].split('.')[0]
                nucs = int(stat_line[4])
                readset_dict[name] = nucs
                all_sample_list.append(name)
    
    ## Get minimum read set            
    min_read_set = min(readset_dict, key=readset_dict.get)
    max_read_set = max(readset_dict, key=readset_dict.get)
    
    ## Restrict minimum
    min_read_set = adjust_min_max(readset_dict, restrict,min_frac)
    min_read_set_size = readset_dict[min_read_set]

    ## Print which read set is smallest
    print('{} was found to be the smallest input read set with {} sequenced nucleotides'.format(min_read_set, min_read_set_size))

    ## Write output to list and string for donwstream processing
    samples = [read_set for read_set in readset_dict.keys() if read_set != min_read_set]
    removed_samples = [read_set for read_set in all_sample_list if read_set not in samples]

    ## Print out which nucleotides are not being used
    #print('The following read sets are smaller than the minimum number of nucleotides for downsampling and will not be processed: {}'.format(removed_samples))

    downsample_nucs = min_read_set_size
    return readset_dict, samples, removed_samples, downsample_nucs

def create_combination_downsamples(readset_dict, read_sets, outlier, minimum):
    number_combos = len(read_sets)

    ## Cases than cannot be handled
    if number_combos > 5:
        sys.exit("You can at most combine five read sets, otherwise the number of possible combinations gets too large")

    if number_combos < 2:
        sys.exit("There are fewer than one read sets available, please inset at least two read sets to be combined/downsampled")

    new_readset_dict = {key: readset_dict[key] for key in read_sets}
    new_min_read_set = min(new_readset_dict, key=new_readset_dict.get)

    if minimum:
        new_min_read_set_size = minimum

    elif outlier and new_min_read_set == outlier:
        tmp_readset_dict = {key: readset_dict[key] for key in read_sets if key != outlier}
        new_min_read_set = min(tmp_readset_dict, key=tmp_readset_dict.get)
        new_min_read_set_size = tmp_readset_dict[new_min_read_set]
    else:
        new_min_read_set_size = new_readset_dict[new_min_read_set]
    
    downsample_nucs = []
    downsample_dict = {}

    ## Cases that can be handled
    if number_combos == 2:
        half_nucs = int(new_min_read_set_size / 2)

        nuc_list = [new_min_read_set_size, half_nucs]
        nuc_dict = {"single" : new_min_read_set_size, "half": half_nucs}
        turned_nuc_dict = {new_min_read_set_size: "single" ,half_nucs: "half"}
        downsample_read_dict = defaultdict(list)
        downsample_read_frac = defaultdict(list)

        for number in nuc_list:
            downsample_nucs.append(number)
            downsample_dict[turned_nuc_dict[number]] = number

        ## downsample operations:
        for read_set in new_readset_dict.keys():
            for number in nuc_list:
                if int(new_readset_dict[read_set]) < int(number):
                    continue
                else:
                    downsample_read_dict[read_set].append(number)
                    downsample_read_frac[read_set].append(turned_nuc_dict[number])

    elif number_combos == 3:
        half_nucs = int(new_min_read_set_size / 2)
        third_nucs = int(new_min_read_set_size / 3)

        nuc_list = [new_min_read_set_size, half_nucs, third_nucs]
        nuc_dict = {"single" : new_min_read_set_size, "half": half_nucs, "third": third_nucs}
        turned_nuc_dict = {new_min_read_set_size: "single" ,half_nucs: "half", third_nucs: "third"}
        downsample_read_dict = defaultdict(list)
        downsample_read_frac = defaultdict(list)

        for number in nuc_list:
            downsample_nucs.append(number)
            downsample_dict[turned_nuc_dict[number]] = number

        ## downsample operations:
        for read_set in new_readset_dict.keys():
            for number in nuc_list:
                if int(new_readset_dict[read_set]) < int(number):
                    continue
                else:
                    downsample_read_dict[read_set].append(number)
                    downsample_read_frac[read_set].append(turned_nuc_dict[number])
    
    elif number_combos == 4:
        half_nucs = int(new_min_read_set_size / 2)
        third_nucs = int(new_min_read_set_size / 3)
        fourth_nucs = int(new_min_read_set_size / 4)

        nuc_list = [new_min_read_set_size, half_nucs, third_nucs, fourth_nucs]
        nuc_dict = {"single" : new_min_read_set_size, "half": half_nucs, "third": third_nucs, "fourth": fourth_nucs}
        turned_nuc_dict = {new_min_read_set_size: "single" ,half_nucs: "half", third_nucs: "third", fourth_nucs: "fourth"}
        downsample_read_dict = defaultdict(list)
        downsample_read_frac = defaultdict(list)

        for number in nuc_list:
            downsample_nucs.append(number)
            downsample_dict[turned_nuc_dict[number]] = number

        ## downsample operations:
        for read_set in new_readset_dict.keys():
            for number in nuc_list:
                if int(new_readset_dict[read_set]) < int(number):
                    continue
                else:
                    downsample_read_dict[read_set].append(number)
                    downsample_read_frac[read_set].append(turned_nuc_dict[number])
    
    elif number_combos == 5:
        half_nucs = int(new_min_read_set_size / 2)
        third_nucs = int(new_min_read_set_size / 3)
        fourth_nucs = int(new_min_read_set_size / 4)
        fifth_nucs = int(new_min_read_set_size / 5)

        nuc_list = [new_min_read_set_size, half_nucs, third_nucs, fourth_nucs, fifth_nucs]
        nuc_dict = {"single" : new_min_read_set_size, "half": half_nucs, "third": third_nucs, "fourth": fourth_nucs, "fifth": fifth_nucs}
        turned_nuc_dict = {new_min_read_set_size: "single" ,half_nucs: "half", third_nucs: "third", fourth_nucs: "fourth", fifth_nucs: "fifth"}
        downsample_read_dict = defaultdict(list)
        downsample_read_frac = defaultdict(list)

        for number in nuc_list:
            downsample_nucs.append(number)
            downsample_dict[turned_nuc_dict[number]] = number

        ## downsample operations:
        for read_set in new_readset_dict.keys():
            for number in nuc_list:
                if int(new_readset_dict[read_set]) < int(number):
                    continue
                else:
                    downsample_read_dict[read_set].append(number)
                    downsample_read_frac[read_set].append(turned_nuc_dict[number])

    return downsample_read_dict, downsample_read_frac, downsample_nucs, downsample_dict


def form_combinations(downsample_dict, downsample_read_frac):
    file_names = [key for key in set(downsample_dict.keys())]
    
    combination_ranges = list(range(1, len(file_names)+1))
    suffixes = [suffix_dict[number] for number in combination_ranges]

    all_combos = {}

    for sf in suffixes:
        #subset_files = [f"{f}.downsampled.{sf}.fastq.gz" for f in file_names]
        subset_files = [f"{f}-{sf}" for f in file_names if sf in downsample_read_frac[f]]
        n = combination_group_sizes[sf]
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

## Not used for now
def main():
    args = parse_args()
    seqkit = args.seqkit
    restrict = args.restrict
    if args.combinations:
        combinations = args.combinations.split(',')
    else:
        print("No combinations are to be done, instead just downsample all reads to the mininum input")

    readset_dict, samples, removed_samples, downsample_nucs = read_seq_stats(seqkit, restrict)
    
    if args.combinations:
        downsample_dict, new_min_read_set, new_min_read_set_size = create_combination_downsamples(readset_dict, combinations)

        form_combinations(downsample_dict)

if __name__ == '__main__':
    main()




