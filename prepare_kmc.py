import argparse
from collections import defaultdict
from itertools import combinations
import sys

DESCRIPTION = ""

def parse_args():
    """Parse CMD args."""

    app = argparse.ArgumentParser(description=DESCRIPTION)
    
    app.add_argument(
    "-u", 
    "--union",
    action="store_true",
    dest="union",
    help="Whether to create a KMC union script. If not invoked, will create difference (or unique kmers if --read_set is invoked)")

    app.add_argument(
    "-i", 
    "--input_kmc",
    action="store",
    nargs='+',
    dest="input",
    help="List of KMC input files, comma separated")

    app.add_argument(
    "-r", 
    "--read_set",
    action="store",
    dest="read_set",
    help="Which read set should be used to compute unique kmers")
    
    app.add_argument(
    "-o", 
    "--output",
    action="store",
    dest="output",
    help="Name of output file (if union or intersection)")
       
    args = app.parse_args()

    return args


def create_union_intersect_script(in_list, out_name, union = True):
    names = []
    with open(out_name, 'w') as out:
        out.write("INPUT:\n")
        for entry in in_list:
            name = entry.split("/")[1].split(".")[0]
            names.append(name)
            out.write("{}={}\n".format(name, "kmc/" + name + '.res'))
        out.write("OUTPUT:\n")
        if union:
            out.write("kmc/kmc_union={}\n".format("+".join(names)))
        else:
            out.write("kmc/kmc_intersection={}\n".format("*".join(names)))

def create_difference_script(read_set, in_list):
    names = []
    out_name_file = "kmc/" + read_set + "_diff_script"
    out_name = "kmc/" + read_set + "_diff"
    with open(out_name_file, 'w') as out:
        out.write("INPUT:\n")
        for entry in in_list:
            name = entry.split("/")[1].split(".")[0]
            if name != read_set:
                names.append(name)
            out.write("{}={}\n".format(name, "kmc/"+ name + '.res'))
        out.write("OUTPUT:\n")
        out.write("{}={} - ({})\n".format(out_name, read_set, "+".join(names)))

def main():
    args = parse_args()
    union = args.union
    list_input = args.input
    read_set = args.read_set
    out_name = args.output

    print(list_input)

    if read_set:
        create_difference_script(read_set, list_input)
    else:
        create_union_intersect_script(list_input, out_name, union)

if __name__ == '__main__':
    main()




