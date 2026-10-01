import re
from collections import defaultdict
import pysam
from tqdm import tqdm 
#import argparse
import pyfastx
import sys

__author__ = "Bernhard Bein, 2024."

# DESCRIPTION = '''
# Script to report read length histograms, partitioned by reads being primary or supplementary alignments
# '''

# def parse_args():
#     """Parse CMD args."""

#     app = argparse.ArgumentParser(description=DESCRIPTION,
#     formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=50, max_help_position=10))
    
#     app.add_argument(
#     "-ib", 
#     "--input_bam",
#     action="store",
#     dest="input_file_bam",
#     help="Name of the input bam file (sorted bam file)")

#     app.add_argument(
#     "-if", 
#     "--input_fastq",
#     action="store",
#     dest="input_file_fastq",
#     help="Name of the input fastq file (fastq.gz)")
    
#     app.add_argument(
#     "-o",
#     "--output_file",
#     action="store",
#     dest="output_file",
#     help="Base name of the output files")
        
#     args = app.parse_args()

#     return args

## This script works directly in snakemake
# def gc_content(seq):
#     """Calculate the GC content of a sequence."""
#     g = seq.count('G')
#     c = seq.count('C')
#     gc_count = g + c
#     return gc_count / len(seq) if len(seq) > 0 else 0

def count_readlen(path_to_fastq):
    fq = pyfastx.Fastx(path_to_fastq)
    readlen_hist = defaultdict(int)
    #gc_hist = defaultdict(lambda: defaultdict(int))
    for name,seq,qual in fq:
        sequence_len = len(seq)
        readlen_hist[sequence_len] += 1
        #gc_content_read = round(gc_content(seq),1)
        #gc_hist[gc_content_read][sequence_len]+=1

    return readlen_hist#,gc_hist

## Helper function to get the absolute length of an alignment through its CIGAR string
def get_len_cigar(cigar_string):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)
    length = 0

    for count, op in matches:
        ## If Cigar string is a deletion relative to the reference
        if op == "D":
            length += 0

        ## if CIGAR operation is not clipping, add length to alignment length
        elif op not in ("S", "H"):
            length += int(count)
    
    return length

## Write a histogram from dictionary
def write_hist_file(readlen_hist, out_path):
    with open(out_path ,'w') as out:
        for key, value in sorted(readlen_hist.items()):
            out.write("{}\t{}\n".format(key, value))

## Write Read GC content
# def write_gc_hist_file(gc_hist, out_path):
#     with open(out_path ,'w') as out:
#         for gc_bin in sorted(gc_hist.keys()):
#             for seq_len in gc_hist[gc_bin].keys():
#                 out.write("{}\t{}\t{}\n".format(gc_bin, seq_len, gc_hist[gc_bin][seq_len]))

## Function to report alignment length, now based on CIGAR-string length
## Might be adapted to use query_alignment_start and query_alignment_end
## Maybe add an option to remove the lenghty "single read output" operation
## Implement mapping types histogram here as well.
def get_aligned_read_length(path_to_bam, outfile):
    bamfile = pysam.AlignmentFile(path_to_bam, "rb")
    total_reads = bamfile.mapped + bamfile.unmapped

    readlen_hist_prim = defaultdict(int)
    readlen_hist_supp = defaultdict(int)
    mapping_types = defaultdict(int)
    with tqdm(total=total_reads, desc="Processing BAM file to get alignment length") as pbar:
        for line in bamfile:
            try:
                read_id = line.query_name
                cigar_string = line.cigarstring
                alignment_type = line.flag
                cigar_len = get_len_cigar(cigar_string)
                ## Add to mapping types histogram:
                mapping_types[alignment_type] += 1
                ## Print alignment length to out
                if alignment_type == 0 or alignment_type == 16:
                    readlen_hist_prim[cigar_len] += 1
                elif alignment_type == 2048 or alignment_type == 2064:
                    readlen_hist_supp[cigar_len] += 1
            except:
                pass

            pbar.update(1)
    out_prim = outfile + "_prim.hist.txt"
    out_sec =  outfile + "_suppl.hist.txt" 
    out_mapping = outfile + "_mappingType.txt" 
    write_hist_file(readlen_hist_prim, out_prim)
    write_hist_file(readlen_hist_supp, out_sec)
    write_hist_file(mapping_types, out_mapping)

def main():
    #args = parse_args()
    infile_bam = snakemake.input[0]#args.input_file_bam
    infile_fastq = snakemake.input[1]#args.input_file_fastq
    outfile = snakemake.params[0]#args.output_file
    #out_hist = snakemake.params[0] + ".hist.txt"
    #out_gc_hist = snakemake.params[0] + ".gc_hist.txt"

    get_aligned_read_length(infile_bam, outfile)

    #hist, gc_hist = count_readlen(snakemake.input[0])
    rl_hist = count_readlen(infile_fastq)
    out_rl_hist = outfile + "_rl.hist.txt"
    write_hist_file(rl_hist, out_rl_hist)   


if __name__ == '__main__':
    main()
