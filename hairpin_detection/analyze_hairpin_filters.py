#!/usr/bin/env python3
import argparse
import itertools
from collections import namedtuple
import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score, recall_score, f1_score, precision_score, roc_auc_score

HairpinInfo = namedtuple("HairpinInfo", ["mapped_bins", "correlation", "slope", "yintercept", "seq_length"])

def pass_mapped_bin_threshold(hairpin_info, args):
    if hairpin_info.mapped_bins < args.mapped_bin_threshold:
        return False
    return True

def pass_correlation(hairpin_info, args):
    if hairpin_info.correlation > args.c:
        return False
    return True

def pass_slope(hairpin_info, args):
    if hairpin_info.slope > args.upper_slope or hairpin_info.slope < args.lower_slope:
        return False
    return True

def pass_yintercept(hairpin_info, args):
    if hairpin_info.yintercept > hairpin_info.seq_length*(1 + args.perc/100) or\
            hairpin_info.yintercept < hairpin_info.seq_length*(1 - args.perc/100):
        return False
    return True

def main():
    parser = argparse.ArgumentParser(description="Detect hairpin artifacts in nanopore reads")
    parser.add_argument("TSV", help="Input ground truth TSV file")
    parser.add_argument("--perc", help="Percentage error allowed for yintercept [10]",
                        type=float, default=10)
    parser.add_argument("--upper_slope", help="Upper threshold for slope [-0.75]", type=float,
                        default=-0.75)
    parser.add_argument("--lower_slope", help="Lower threshold for slope [-1.25]", type=float,
                        default=-1.25)
    parser.add_argument("-c", help="Threshold for correlation [-0.75]", type=float,
                        default=-0.75)
    parser.add_argument("-m", "--mapped-bin-threshold",
                        help="Threshold number of bins with mapped minimizers [5]",
                        type=int, default=5)

    args = parser.parse_args()

    df, filter_combinations_str, filter_keys = run_all_filters(args)

    print("Filter", "\t".join(filter_keys), "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", sep="\t")

    for filter_combination in filter_combinations_str:
        tn, fp, fn, tp = confusion_matrix(df["Hairpin_status"], df[filter_combination],
                                          labels=["Non-hairpin", "Hairpin"]).ravel()
        filter_combo_summary = "\t".join(["X" if filt in filter_combination.split("_") else " " for filt in filter_keys ])
        print(filter_combination, filter_combo_summary,
              accuracy_score(df["Hairpin_status"], df[filter_combination]),
              recall_score(df["Hairpin_status"], df[filter_combination], pos_label="Hairpin"),
              (tn / (tn + fp)),
              precision_score(df["Hairpin_status"], df[filter_combination], pos_label="Hairpin"),
              f1_score(df["Hairpin_status"], df[filter_combination], pos_label="Hairpin"), sep="\t")


def run_all_filters(args):
    filters = {"mappedbins": pass_mapped_bin_threshold, "yint": pass_yintercept,
               "slope": pass_slope, "correlation": pass_correlation}
    filter_keys = list(filters.keys())
    filter_combinations = []
    for i in range(1, len(filters) + 1):
        filter_combinations += itertools.combinations(filters.keys(), i)
    filter_combinations_str = ["_".join(el) for el in filter_combinations]
    print(filter_combinations_str)
    df_content = []
    with open(args.TSV, 'r') as fin:
        fin.readline()
        for line in fin:
            ctg_line = []
            line = line.strip().split("\t")
            name, hairpin_ground_truth = line[:2]
            length, correlation, yint, slope = line[4:8]
            mapped_bins, hairpin_pred = line[9:11]

            hairpin_info = HairpinInfo(mapped_bins=int(mapped_bins), correlation=float(correlation),
                                       slope=float(slope), yintercept=float(yint), seq_length=int(length))

            ctg_line.extend([name, hairpin_ground_truth])

            for filter_combination in filter_combinations:
                hairpin_pred = "Hairpin"
                for hp_filter in filter_combination:
                    if not filters[hp_filter](hairpin_info, args):
                        hairpin_pred = "Non-hairpin"
                        break
                ctg_line.append(hairpin_pred)
            df_content.append(ctg_line)
    columns_list = ["Name", "Hairpin_status"] + filter_combinations_str
    df = pd.DataFrame(df_content, columns=columns_list)
    return df, filter_combinations_str, filter_keys


if __name__ == "__main__":
    main()
