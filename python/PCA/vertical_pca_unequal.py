"""
    Copyright (C) 2020 Anne Hartebrodt

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    Authors: Anne Hartebrodt

"""

# import import_export.easy_import as easy
import argparse as ap
import os
import os.path as op
import time
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.sparse.linalg as lsa
from scipy.sparse import coo_matrix
import sys

import python.PCA.comparison  as co
import python.PCA.convenience as cv
import python.PCA.shared_functions as sh
import python.PCA.vertical_pca_library as gv
import python.PCA.vertical_pca_runner as runner
import python.import_export.gwas_import as gi
import python.import_export.mnist_import as mi
import python.import_export.spreadsheet_import as si

from python.PCA.vertical_pca_benchmark import *

def unequal_benchmark_edge_case(nr_splits=3):
    splits = []
    for p in range(1, 30, 1):
        p1 = p * 0.01
        p2 = (1.0 - p * 0.01) / 2
        splits.append([p1, p2, p2])
    splits = np.asarray(splits)
    return splits

def run_unequal_benchmark(data, dataset_name, maxit, nr_repeats, k,outdir, splits=[[0.01, 0.99]], convergence_eps=1e-6, precomputed_pca=None, nr_samples=-1, nr_features=-1):
    '''
    run the simulation of a federated run of vertical power iteration
    Args:
        data: data frame or list of data frames containing dimension which is split
        in the columns
        dataset_name: Name, for logging
        maxit: maximal iterations to run
        nr_repeats: number of times to repeat experiments
        k: targeted dimensions
        splits: array of splits for dataset (only applicable when data is not list)
        outdir: result directory

    Returns:

    '''
    # g = gv.standalone(data, k)
    islist = False
    if isinstance(data, list):
        islist = True
        data_list = data
        data = np.concatenate(data, axis=1)
        splits = [1]  # data is already split, only counter experiments need to be run.

    u, s, v = lsa.svds(data.T, k=k)
    u = np.flip(u, axis=1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    dataset_name_guo = dataset_name + '_guo'

    for c in range(nr_repeats):
        for s in splits:
            # filename will be the same for angle log file and correlation log file
            timer = time.monotonic()


            # split the data
            if not islist:
                data_list, choice = sh.partition_data_vertically(data, s, randomize=True, equal=False, splits=s)
            else:
                choice = range(data.shape[1])
            logf = op.join(outdir, 'log_choices.log')
            log_choices(logf, filename, choice)

            start = time.monotonic()
           #logftime = op.join(outdir, 'time.log')

            # # simulate the run
            for fedqr, mode in zip([True, False], ['fed_qr', 'central_qr']):
                filename = init_benchmark(outdir=outdir, dataset_name=dataset_name+'_'+mode, maxit=maxit, counter=c,
                                          nr_samples=nr_samples, nr_features=nr_features, k=k,
                                          convergence_eps=convergence_eps,
                                          splits=s, timer=timer, transmission_costs=True)
                simulate_guo_benchmark(data_list, k, maxit=maxit, scipy=u, filename=filename, choices=choice,
                                       precomputed_pca=precomputed_pca, federated_qr=fedqr)
                end = time.monotonic()

               #log_time(logftime, 'qr_scheme'+'_'+mode, end - start, s, c)

                filename = init_benchmark(outdir=outdir, dataset_name=dataset_name_guo+'_'+mode, maxit=maxit, counter=c,
                                          nr_samples=nr_samples, nr_features=nr_features, k=k,
                                          convergence_eps=convergence_eps,
                                          splits=s, timer=timer, transmission_costs=True)

                start = time.monotonic()
                compute_k_eigenvectors(data_list, k=k, maxit=maxit, scipy=u, filename=filename, choices=choice,
                                                precomputed_pca=precomputed_pca, federated_qr=fedqr)
                end = time.monotonic()
               #log_time(logftime, 'guo_single'+'_'+mode, end - start, s, c)

                filename = init_benchmark(outdir=outdir, dataset_name=dataset_name + '_hybrid'+'_'+mode, maxit=maxit, counter=c,
                                          nr_samples=nr_samples, nr_features=nr_features, k=k,
                                          convergence_eps=convergence_eps,
                                          splits=s, timer=timer, transmission_costs=True)
                filename2 = init_benchmark(outdir=outdir, dataset_name=dataset_name + '_hybrid_reit'+'_'+mode, maxit=maxit,
                                           counter=c,
                                           nr_samples=nr_samples, nr_features=nr_features, k=k,
                                           convergence_eps=convergence_eps,
                                           splits=s, timer=timer, transmission_costs=True)

                start = time.monotonic()
                hybrid_scheme(data_list, k=k, maxit=maxit, scipy=u, filename=filename, filename2=filename2,
                             choices=choice,
                             precomputed_pca=precomputed_pca, federated_qr=fedqr)
                end = time.monotonic()
               #log_time(logftime, 'hybrid_scheme'+'_'+mode, end - start, s, c)


if __name__ == '__main__':
    #data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')

    data, sample_ids, variable_names = si.data_import('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t', header=0, rownames=0)
    data = si.scale_center_data_columnwise(data)
    #data = coo_matrix.asfptype(data)
    # args.k = 10
    # g = gv.standalone(data, k=2)
    #
    u, s, v = lsa.svds(data.T,k=10)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    splits =    [[0.01 , 0.495, 0.495],
                [0.01, 0.01, 0.98],
                [0.025, 0.025, 0.05, 0.9]]

