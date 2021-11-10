#!/usr/bin/env python3
'''
Computing statistics on minimizers  from reads
'''
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

def compute_read_statistics(mx_info):
    "Compute various statistics on the given minmizer sketch of the read"
    mx_df = make_dataframe(mx_info)
    pearson_corr = pearson_correlation_coefficient(mx_df)
    spearman_corr = spearman_correlation_coefficient(mx_df)
    yintercept, slope = robust_linear_regression(mx_df)

    return pearson_corr, spearman_corr, yintercept, slope
