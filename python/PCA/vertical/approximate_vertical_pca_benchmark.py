from python.PCA.vertical.vertical_pca_benchmark import simulate_subspace_iteration

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
import python.import_export.gwas_import as gi
import python.import_export.mnist_import as mi
import python.import_export.spreadsheet_import as si
import python.PCA.vertical.federated_qr as qr
import python.PCA.horizontal.power_iteration as powerit
import json

from python.PCA.horizontal.horizontal_pca_power_iteration import simulate_distributed_horizontal
from python.PCA.logging import *
from python.PCA.horizontal.balcan import simulate_federated_horizontal_pca

## Approximate federated vertical PCA
## Speed up over regular version
def benchmark_vertical_approximate_pca(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-9, unequal=False, precomputed_pca=None):
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

    current_split = 0
    for c in range(nr_repeats):
        for s in splits:
            # split the data
            if not islist:
                if unequal:
                    # if unequal then the number of sites is the length of s and s itself contains the splits
                    data_list, choice = sh.partition_data_vertically(data, len(s), randomize=True, perc=s, equal=False)
                    s = current_split
                    current_split += 1
                else:
                    data_list, choice = sh.partition_data_vertically(data, s, randomize=True)
            else:
                # make a dummy choice vector
                choice = range(data.shape[1])

            logftime = op.join(outdir, 'time.log')

            # # simulate the run
            start = time.monotonic()
            # simultaneous only H
            grad = False
            grad_name = 'power'
            mode = 'complete'
            ortho_freq = 1000 # extremely high, will rarely be reached
            fedqr = False
            print('power - matrix - ' + mode)
            outdir_gradient = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_gradient, exist_ok=True)
            filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)



            # Compute federated approximate PCA
            start = time.monotonic()
            mode = 'approximative'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            g, h = approximate_vertical(data_list, k, factor_k=2)
            g = np.concatenate(g, axis=0)
            log_current_accuracy(u=u, G_i=g, current_iteration=1,
                                 filename=filename, precomputed_pca=precomputed_pca, v=v, H_i=h, choices= choice)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)

            # Compute full decomposition using approximative PCA as seed
            start = time.monotonic()
            print('test')
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k])
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)


def approximate_vertical(data_list, k=10, factor_k=2):
    data_list = [d.T for d in data_list]
    v, e = simulate_federated_horizontal_pca(data_list, k=k, factor_k=factor_k)
    g = [np.dot(d, v) for d in data_list]
    return g, v



if __name__ == '__main__':
    print('test')

    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)

    dataset_name = 'mnist'
    maxit = 500
    nr_repeats = 10
    k = 10
    splits = [5, 10]
    outdir = '/home/anne/Documents/featurecloud/pca/approximative-vertical/results'
    benchmark_vertical_approximate_pca(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-9,
                                           unequal=False, precomputed_pca=None)
