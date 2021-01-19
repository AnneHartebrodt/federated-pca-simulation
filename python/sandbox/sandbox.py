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
import python.PCA.federated_qr as fqr

####### MATRIX POWER ITERATION SCHEME #######
def simulate_guo_benchmark(local_data, k, maxit, filename=None, scipy=None, choices=None, precomputed_pca=None, fractev=1.0, encrypt=False):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = []
    iterations = 0
    converged = False
    total_len = 0
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    # send parts to local sites
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        #log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
        start = start + local_data[i].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    converged_eigenvals = []
    while not converged and iterations < maxit and len(converged_eigenvals)<k*fractev:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], k))
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            #log_transmission(filename, "H_local=CS", iterations, i, H_local)
            H_i = H_i + H_local
        #log_transmission(filename, "H_global=SC", iterations, 1, H_i)

        for i in range(len(G_list)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            #log_transmission(filename, "Gi_local=CS", iterations, i, G_list[i])

        G_i = np.concatenate(G_list, axis=0)

        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        eigenvals = np.sqrt(eigenvals)

        G_i, R = la.qr(G_i, mode='economic')

        orho, G_list = fqr.simulate_federated_qr(G_list, encrypt=encrypt)
        # start = 0
        # for i in range(len(G_list)):
        #     G_list[i] = G_i[start:start + local_data[i].shape[1], :]
        #     #log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
        #     start = start + local_data[i].shape[1]
        print(co.compute_angles(orho, orho[:,1:]))
        converged, conv, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)
        H_i_prev = H_i
        #log_current_accuracy(scipy=scipy, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations,
               #              filename=filename,
                    #         choices=choices, precomputed_pca=precomputed_pca)
    G_i = np.concatenate(G_list)
    print(iterations)
    return G_i, eigenvals, converged_eigenvals

if __name__ == '__main__':
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')

    #data, sample_ids, variable_names = si.data_import('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t', header=0, rownames=0)
    #data = si.scale_center_data_columnwise(data)
    data = coo_matrix.asfptype(data)
    # args.k = 10
    # g = gv.standalone(data, k=2)
    #
    u, s, v = lsa.svds(data.T,k=10)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    data_list, choices = sh.partition_data_vertically(data,perc=[10,30,60], equal=False)

    #ev =compute_k_eigenvectors(data_list, 10, 2000)
    #ev1 = hybrid_scheme(data_list, 10, 2000, scipy=u)
    start = time.monotonic()
    #ev2,ee2,ee3 = simulate_guo_benchmark(data_list, 7, 500, encrypt=False)
    end = time.monotonic()
    print(end-start)

    start = time.monotonic()
    ev2,ee2,ee3 = simulate_guo_benchmark(data_list, 7, 500, encrypt=False)
    end = time.monotonic()
    print(end-start)
    #print(co.compute_angles(u, ev))
    #print(co.compute_angles(u, ev1))
    print(co.compute_angles(u, ev2))
   # print(co.angle(u[:,7], ev[:,9]))