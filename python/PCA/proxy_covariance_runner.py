import argparse as ap
import pandas as pd
import import_export.easy_import as easy
import proxy_covariance as pca
import convenience as cv
import balcan as bpca
import numpy as np

def simulateImtiaz(datasets, totaln):
    partial = []
    for d in datasets:
        # same function for balcan and imtiaz
        partial.append(bpca.perform_SVD(d))
    weights = [d.shape[0]/totaln for d in datasets]
    dpca = pca.aggregate_partial_SVDs(partial, 6, weights)
    return dpca

def simulateQu(datasets, totaln):
    partial = []
    sums = []
    for d in datasets:
        # same function for balcan and imtiaz
        partial.append(bpca.perform_SVD(d))
        sums.append(np.sum(d, axis=0))
    weights = [d.shape[0]/totaln for d in datasets]
    dpca = pca.aggregate_partial_SVDs(partial, 6, weights, sums, totaln)
    return dpca

