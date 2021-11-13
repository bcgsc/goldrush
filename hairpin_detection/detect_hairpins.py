#!/usr/bin/env python3
'''
Detect hairpin structure artifacts in oxford nanopore reads
'''
import argparse
from collections import defaultdict, namedtuple, Counter
import igraph as ig
import ntlink_utils
import sys
import calculate_hairpin_stats

MinimizerInfo = namedtuple("MinimizerInfo", ["pos", "strand", "time_seen"])


def tally_sequence_lengths(index_filename):
    "Tally the sequence lengths from the faidx index file"
    lengths = {}
    with open(index_filename, 'r') as fin:
        for line in fin:
            name, length = line.strip().split("\t")[:2]
            lengths[name] = int(length)
    return lengths

def is_valid_mx(list_mx_info):
    "Return true if the given minimizer has a multiplicity of 2 and those are found from different strands"
    if len(list_mx_info) != 2:
        return False
    info1, info2 = list_mx_info
    if info1.strand == info2.strand:
        return False
    return True

def is_valid_position(pos, end_length, seq_length):
    "Return True if position falls within prescribed end lengths"
    if pos <= end_length or pos >= (seq_length - end_length):
        return True
    return False

def has_valid_positions(minimizer_info, args, length):
    "Return True if both instances of minimizer have valid positions"
    assert len(minimizer_info) == 2
    return is_valid_position(minimizer_info[0].pos, args.e, length) and \
           is_valid_position(minimizer_info[1].pos, args.e, length)


def filter_ordered_sketch(mx_line, args, seq_length):
    "Given a line of a indexlr file, parse and only keep minimizers with multiplicity of 2, strictly on different strands"

    mx_info = defaultdict()  # mx -> [Minimizer info]
    line = mx_line.strip().split("\t")
    if len(line) > 1:
        name, mxs_all = line
        mx_pos_strands = mxs_all.split(" ")
        for mx_pos_strand in mx_pos_strands:
            mx, pos, strand = mx_pos_strand.split(":")
            # if not is_valid_position(int(pos), args.e, seq_length):
            #     continue
            if mx not in mx_info:
                mx_info[mx] = [MinimizerInfo(int(pos), strand, 1)] #!! TODO: change numbers?
            else:
                mx_info[mx].append(MinimizerInfo(int(pos), strand, 2))

        mx_info = {mx: mx_info[mx] for mx in mx_info if is_valid_mx(mx_info[mx]) and
                   has_valid_positions(mx_info[mx], args, seq_length)}
        mxs = [mx_pos_strand.split(":")[0] for mx_pos_strand in mx_pos_strands
               if mx_pos_strand.split(":")[0] in mx_info]
        for mx in mx_info:
            assert len(mx_info[mx]) == 2
        assert len(mxs) == len(mx_info.keys())*2

    return mx_info, mxs


def is_hairpin(mx_info, correlation, yintercept, slope, seq_length, args):
    "Return true if fits the threshold for a hairpin"
    if len(mx_info) < 3:
        return False
    if correlation > args.c:
        return False
    if slope > args.upper_slope or slope < args.lower_slope:
        return False
    if yintercept > seq_length*1.10 or yintercept < seq_length*0.9:
        return False
    return True

def detect_hairpins(args, seq_lengths):
    "Read through minimizers for each read, and determine whether it is a putative hairpin artifact"
    hairpins = 0
    total_reads = 0

    fout = open(args.o, 'w')
    fout.write("Name\tLength\tCorrelation_coefficient\tyintercept\tslope\tnum_mx\tis_hairpin_pred\n")

    format_str = ("{}\t"*7).strip() + "\n"

    with open(args.MX, 'r') as mx_in:
        for mx_line in mx_in:
            name, _ = mx_line.strip().split("\t")
            mx_info, mxs = filter_ordered_sketch(mx_line, args, seq_lengths[name])

            correlation, yint, slope = 0, 0, 0
            if len(mx_info) >= 3:
                correlation, yint, slope = calculate_hairpin_stats.compute_read_statistics(mx_info, args.corr)

            if is_hairpin(mx_info, correlation, yint, slope, seq_lengths[name], args):
                hairpins += 1
                fout.write(format_str.format(name, seq_lengths[name], correlation, yint, slope, len(mxs), "Hairpin"))
            else:
                fout.write(format_str.format(name, seq_lengths[name], correlation, yint, slope, len(mxs), "Non-hairpin"))

            total_reads += 1

    fout.close()
    return hairpins, total_reads


def print_args(args):
    "Print the values of the arguments"
    print("\nHairpin detection parameters:")
    print("\tMX", args.MX)
    print("\t-i", args.index)
    print("\t--perc", args.perc)
    print("\t-e", args.e)
    print("\t--upper_slope", args.upper_slope)
    print("\t--lower_slope", args.lower_slope)
    print("\t-c", args.c)
    print("\t--corr", args.corr)


def main():
    "Detect hairpin structures in input nanopore reads from minimizer sketches"
    parser = argparse.ArgumentParser(description="Detect hairpin artifacts in nanopore reads")
    parser.add_argument("MX", help="Input minimizers TSV file, or '-' if piping to standard in")
    parser.add_argument("-i", "--index", help="samtools faidx index for input reads", required=True, type=str)
    parser.add_argument("--perc", help="Percentage error allowed for yintercept", type=float, default=10)
    parser.add_argument("-e", help="Length of ends to consider", type=int, default=5000)
    parser.add_argument("--upper_slope", help="Upper threshold for slope", type=float, default=-0.75)
    parser.add_argument("--lower_slope", help="Lower threshold for slope", type=float, default=-1.25)
    parser.add_argument("-c", help="Threshold for correlation", type=float, default=-0.75)
    parser.add_argument("--corr", help="Correlation coefficient to use. Valid values are pearson or spearman",
                        default="spearman", type=str)
    parser.add_argument("-o", help="Output file for hairpin classifications [stdout]", type=str, default=sys.stdout)
    args = parser.parse_args()

    print("Running hairpin detection...")

    print_args(args)

    if args.corr not in ["pearson", "spearman"]:
        raise ValueError("--corr must be set to pearson or spearman. ", args.corr, "supplied.")

    args.MX = "/dev/stdin" if args.MX == "-" else args.MX

    seq_lengths = tally_sequence_lengths(args.index)
    hairpins, total_reads = detect_hairpins(args, seq_lengths)
    print("Total reads analyzed:", total_reads)
    print("Total hairpins detected:", hairpins)

    print("DONE!")


if __name__ == "__main__":
    main()
