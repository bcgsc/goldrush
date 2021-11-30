#!/usr/bin/env python3
'''
Computing statistics on minimizers  from reads
'''
from collections import Counter
from scipy.stats import entropy
import statsmodels.api as sm
import pandas as pd
import scipy.stats as sci


def make_dataframe(mx_info):
    "Given a dictionary of mx -> [MinimizerInfo], generate a pandas dataframe"
    formatted_data = [(mx_info[mx][0].pos, mx_info[mx][1].pos) for mx in mx_info]
    df = pd.DataFrame(formatted_data, columns=["position1", "position2"])
    return df

def pearson_correlation_coefficient(df):
    "Return Pearson correlation coefficient of dataframe columns"
    return sci.pearsonr(df["position1"], df["position2"])[0]

def spearman_correlation_coefficient(df):
    "Return Spearman correlation coefficient of dataframe columns"
    return sci.spearmanr(df["position1"], df["position2"]).correlation

def robust_linear_regression(df):
    "Perform robust linear regression on dataframe columns"
    X = df[["position1"]]
    y = df[["position2"]]

    rlm_model = sm.RLM(y, sm.add_constant(X))
    rlm_results = rlm_model.fit()

    p = rlm_results.params

    return p.const, p.position1

def find_correlation_coefficient(mx_df, correlation_arg):
    if correlation_arg == "pearson":
        return pearson_correlation_coefficient(mx_df)
    else:
        return spearman_correlation_coefficient(mx_df)

def compute_entropy(mx_df, read_len, end_len):
    end_len = int(read_len/2) if 2*end_len > read_len else end_len
    bins = pd.cut(range(0, end_len+1), bins=10, retbins=True)[1]
    max_entropy = entropy([1/10]*10, base=2)
    cut_bins = pd.cut(mx_df["position1"], bins=bins)
    counts_bins = Counter(cut_bins[0]) # TODO: magic number
    cut_bins_bins = counts_bins[0].categories
    for bin_int in cut_bins_bins:
        if bin_int not in counts_bins:
            counts_bins[bin_int] = 0

    all_rows = len(mx_df)
    entropy_hairpin = entropy([counts_bins[int_bin]/all_rows for int_bin in counts_bins], base=2)

    chi_test = sci.chisquare([counts_bins[int_bin] for int_bin in counts_bins]).pvalue
    chi_test_expected = sci.chisquare([counts_bins[int_bin] for int_bin in counts_bins], [int(all_rows/10)]*10).pvalue

    return entropy_hairpin/max_entropy, chi_test, chi_test_expected # TODO How to have empty bins? Divide the read in half?

def compute_read_statistics(mx_info, args, read_len):
    "Compute various statistics on the given minimizer sketch of the read"
    mx_df = make_dataframe(mx_info)
    corr = find_correlation_coefficient(mx_df, args.corr)
    yintercept, slope = robust_linear_regression(mx_df)
    entropy_calc, chi1, chi2 = compute_entropy(mx_df, read_len, args.e)

    return corr, yintercept, slope, entropy_calc, chi1, chi2
