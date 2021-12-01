#!/usr/bin/env python3
'''
Computing statistics on minimizers  from reads
'''
from collections import Counter
import warnings
from scipy.stats import entropy
import statsmodels.api as sm
import statsmodels.tools.sm_exceptions as sm_except
import pandas as pd
import scipy.stats as sci
warnings.simplefilter(action='ignore', category=sm_except.ConvergenceWarning)


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
    "Returns specified correlation coefficient"
    if correlation_arg == "pearson":
        return pearson_correlation_coefficient(mx_df)
    return spearman_correlation_coefficient(mx_df)

def compute_entropy(mx_df, read_len, end_len, num_bins=10):
    "Compute entropy stats"
    end_len = int(read_len/2) if 2*end_len > read_len else end_len
    bins = pd.cut(range(0, end_len+1), bins=num_bins, retbins=True)[1]
    max_entropy = entropy([1/num_bins]*num_bins, base=2)
    cut_bins = pd.cut(mx_df["position1"], bins=bins, retbins=True)
    counts_bins = Counter(cut_bins[0])
    cut_bins_bins = cut_bins[0].cat.categories
    num_mapped_bins = len(counts_bins)
    for bin_int in cut_bins_bins:
        if bin_int not in counts_bins:
            counts_bins[bin_int] = 0

    all_rows = len(mx_df)
    entropy_hairpin = entropy([counts_bins[int_bin]/all_rows for int_bin in counts_bins], base=2)

    return entropy_hairpin/max_entropy, num_mapped_bins

def compute_read_statistics(mx_info, args, read_len):
    "Compute various statistics on the given minimizer sketch of the read"
    mx_df = make_dataframe(mx_info)
    corr = find_correlation_coefficient(mx_df, args.corr)
    yintercept, slope = robust_linear_regression(mx_df)
    entropy_calc, mapped_bins = compute_entropy(mx_df, read_len, args.e, num_bins=args.bins)

    return corr, yintercept, slope, entropy_calc, mapped_bins
