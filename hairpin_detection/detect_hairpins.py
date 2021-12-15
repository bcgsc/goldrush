#!/usr/bin/env python3
'''
Detect hairpin structure artifacts in oxford nanopore reads
'''
import argparse
from collections import defaultdict, namedtuple
import sys
from joblib import load
import btllib
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
    "Return true if the given minimizer has a multiplicity of 2, from different strands"
    if len(list_mx_info) != 2:
        return False
    info1, info2 = list_mx_info
    if info1.strand == info2.strand:
        return False
    return True


def is_valid_position(mx_info, args, seq_length):
    "Return True if position falls within prescribed end lengths"
    if 2*args.e > seq_length:
        end_length = int(seq_length/2)
    else:
        end_length = args.e
    for mi in mx_info:
        if mi.time_seen == 1:
            if mi.pos > end_length:
                return False
        elif mi.time_seen == 2:
            if mi.pos < (seq_length - end_length):
                return False
        else:
            raise ValueError("time_seen should be 1 or 2")
    return True


def has_valid_positions(minimizer_info, args, seq_length):
    "Return True if both instances of minimizer have valid positions"
    assert len(minimizer_info) == 2
    return is_valid_position(minimizer_info, args, seq_length)


def filter_ordered_sketch(mxs, args, seq_length):
    "Given minimizer sketch, parse and keep minimizers with multiplicity of 2 on different strands"

    mx_info = defaultdict()  # mx -> [Minimizer info]
    mx_track = defaultdict() # mx -> MinimizerInfo
    if len(mxs) > 1:
        for mx_pos_strand in mxs:
            mx, pos, strand = mx_pos_strand.out_hash, mx_pos_strand.pos, mx_pos_strand.forward
            if mx not in mx_track:
                mx_track[mx] = MinimizerInfo(int(pos), strand, 1)
            else:
                if mx not in mx_info:
                    mx_info[mx] = [mx_track[mx]]
                mx_info[mx].append(MinimizerInfo(int(pos), strand, 2))

        mx_info = {mx: mx_info[mx] for mx in mx_info if is_valid_mx(mx_info[mx]) and
                   has_valid_positions(mx_info[mx], args, seq_length)}
        for mx in mx_info:
            assert len(mx_info[mx]) == 2

    return mx_info


def is_hairpin(mx_info, correlation, yintercept, slope, mapped_bins, seq_length, args):
    "Return true if fits the threshold for a hairpin"
    if len(mx_info) < args.mapped_bin_threshold:
        return False
    if correlation > args.c:
        return False
    if slope > args.upper_slope or slope < args.lower_slope:
        return False
    if yintercept > seq_length*(1 + (args.perc/100)) or yintercept < seq_length*(1 - (args.perc/100)):
        return False
    if mapped_bins < args.mapped_bin_threshold:
        return False
    return True

def detect_hairpins(args, seq_lengths):
    "Read through minimizers for each read, and determine whether it is a putative hairpin artifact"
    hairpins = 0
    total_reads = 0

    fout = open(args.o, 'w')
    fout.write("Name\tLength\tCorrelation_coefficient\tyintercept\tslope\tnum_mx"
               "\tmapped_bins\tis_hairpin_pred\trf_pred\n")

    format_str = ("{}\t"*9).strip() + "\n"

    # Load models for random forest
    if args.r:
        classifier = load(args.r + "/random_forest_classifier")
        scaler = load(args.r + "/scaler")

    with btllib.Indexlr(args.FA, args.k, args.w, btllib.IndexlrFlag.LONG_MODE, 2) as minimizers: # !! TODO: specify flags when bug fixed
        for mx_entry in minimizers:
            name = mx_entry.id
            mx_info = filter_ordered_sketch(mx_entry.minimizers, args, seq_lengths[name])

            if args.v:
                print("Name", "Minimizer1", "Minimizer2", sep="\t", file=sys.stderr)
                for mx in mx_info:
                    assert len(mx_info[mx]) == 2
                    mx_list = mx_info[mx]
                    assert mx_list[0].pos < mx_list[1].pos
                    print(name, mx_list[0].pos, mx_list[1].pos, sep="\t", file=sys.stderr)

            correlation, yint, slope, mapped_bins = 0, 0, 0, 0
            if args.r:
                random_forest_classification = "Non-hairpin"
            else:
                random_forest_classification = "N/A"

            if len(mx_info) >= args.mapped_bin_threshold:
                correlation, yint, slope, mapped_bins = \
                    calculate_hairpin_stats.compute_read_statistics(mx_info, args,
                                                                    seq_lengths[name])
                if args.r:
                    random_forest_classification = \
                        calculate_hairpin_stats.random_forest(correlation, slope, len(mx_info),
                                                              mapped_bins, seq_lengths[name]/yint,
                                                              classifier, scaler)

            if is_hairpin(mx_info, correlation, yint, slope, mapped_bins, seq_lengths[name], args):
                hairpins += 1
                fout.write(format_str.format(name, seq_lengths[name], correlation, yint, slope,
                                             len(mx_info), mapped_bins, "Hairpin",
                                             random_forest_classification))
            else:
                fout.write(format_str.format(name, seq_lengths[name], correlation, yint, slope,
                                             len(mx_info), mapped_bins, "Non-hairpin",
                                             random_forest_classification))

            total_reads += 1

    fout.close()
    return hairpins, total_reads


def print_args(args):
    "Print the values of the arguments"
    print("\nHairpin detection parameters:")
    print("\tFA", args.FA)
    print("\t-i", args.index)
    print("\t--perc", args.perc)
    print("\t-e", args.e)
    print("\t--upper_slope", args.upper_slope)
    print("\t--lower_slope", args.lower_slope)
    print("\t-c", args.c)
    print("\t-b", args.bins)
    print("\t-m", args.mapped_bin_threshold)
    print("\t--corr", args.corr)


def main():
    "Detect hairpin structures in input nanopore reads from minimizer sketches"
    parser = argparse.ArgumentParser(description="Detect hairpin artifacts in nanopore reads")
    parser.add_argument("FA", help="Input fasta file, or '-' if piping to standard in")
    parser.add_argument("-i", "--index", help="samtools faidx index for input reads",
                        required=True, type=str)
    parser.add_argument("-k", help="Kmer size", required=True, type=int)
    parser.add_argument("-w", help="Window size", required=True, type=int)
    parser.add_argument("--perc", help="Percentage error allowed for yintercept [10]",
                        type=float, default=10)
    parser.add_argument("-e", help="Length of ends to consider (bp) [5000]", type=int, default=5000)
    parser.add_argument("--upper_slope", help="Upper threshold for slope [-0.75]", type=float,
                        default=-0.75)
    parser.add_argument("--lower_slope", help="Lower threshold for slope [-1.25]", type=float,
                        default=-1.25)
    parser.add_argument("-c", help="Threshold for correlation [-0.75]", type=float,
                        default=-0.75)
    parser.add_argument("--corr", help="Correlation coefficient to use. "
                                       "Valid values are pearson or spearman [spearman]",
                        default="spearman", type=str)
    parser.add_argument("-b", "--bins", help="Number of bins for minimizer distribution check [10]",
                        type=int, default=10)
    parser.add_argument("-m", "--mapped-bin-threshold",
                        help="Threshold number of bins with mapped minimizers [5]",
                        type=int, default=5)
    parser.add_argument("-o", help="Output file for hairpin classifications [stdout]",
                        type=str, default=sys.stdout)
    parser.add_argument("-r", help="Path to random forest models", required=False)
    parser.add_argument("-v", action="store_true", help="Verbose logging of filtered minimizers")

    args = parser.parse_args()

    print("Running hairpin detection...")

    print_args(args)

    if args.corr not in ["pearson", "spearman"]:
        raise ValueError("--corr must be set to pearson or spearman. ", args.corr, "supplied.")

    args.FA = "/dev/stdin" if args.FA == "-" else args.FA

    seq_lengths = tally_sequence_lengths(args.index)
    hairpins, total_reads = detect_hairpins(args, seq_lengths)
    print("Total reads analyzed:", total_reads)
    print("Total hairpins detected:", hairpins)

    print("DONE!")


if __name__ == "__main__":
    main()
