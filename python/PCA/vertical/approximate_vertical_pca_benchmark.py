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
from python.PCA.horizontal.balcan import simulate_federated_horizontal_pca, local_SVD

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
            start = time.monotonic()
            # simultaneous only H
            grad = False
            grad_name = 'power'
            mode = 'complete'
            ortho_freq = maxit+1 # will not be reached
            fedqr = False

            # # # simulate the run
            start = time.monotonic()
            print('test')
            mode = 'randomized-1'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)

            run_randomized(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
                           precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                           epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))
            start = time.monotonic()
            print('test')
            mode = 'randomized-2'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)

            run_randomized(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
                           precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                           epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))
            start = time.monotonic()
            print('test')
            mode = 'randomized-3'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)

            run_randomized(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
                           precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                           epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None, use_approximate=False)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))

            print('power - matrix - ' + mode)
            outdir_gradient = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_gradient, exist_ok=True)
            filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))


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
            print(mode + ' ' + str(end - start))


            #Compute full decomposition using approximative PCA as seed
            start = time.monotonic()
            print('test')
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k])
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            #
            print(mode + ' ' + str(end - start))
            # Compute federated approximate PCA
            start = time.monotonic()
            mode = 'approximative-smpc-add'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            g, h= add_approx_vertical(data_list, k, factor_k=2)
            g = np.concatenate(g, axis=0)
            log_current_accuracy(u=u, G_i=g, current_iteration=1,
                                 filename=filename, precomputed_pca=precomputed_pca, v=v, H_i=h, choices=choice)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))
            #Compute full decomposition using approximative PCA as seed
            start = time.monotonic()
            print('test')
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k])
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))
            # Compute federated approximate PCA
            start = time.monotonic()
            mode = 'approximative-smpc'
            outdir_approx = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            g, h = approximate_vertical_smpc(data_list, k, factor_k=2)
            g = np.concatenate(g, axis=0)
            log_current_accuracy(u=u, G_i=g, current_iteration=1,
                                 filename=filename, precomputed_pca=precomputed_pca, v=v, H_i=h, choices=choice)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode+' '+str(end-start))

            # Compute full decomposition using approximative PCA as seed
            start = time.monotonic()
            print('test')
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k])
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)
            print(mode + ' ' + str(end - start))

            #


def approximate_vertical(data_list, k=10, factor_k=2):
    data_list = [d.T for d in data_list]
    v, e = simulate_federated_horizontal_pca(data_list, k=k, factor_k=factor_k)
    g = [np.dot(d, v) for d in data_list]
    return g, v

def approximate_vertical_smpc(data_list, k=10, factor_k=2):
    data_list = [d.T for d in data_list]
    v = smpc_enabled_approx_vertical(data_list, k=k, factor_k=factor_k)
    g = [np.dot(d, v) for d in data_list]
    return g, v

def smpc_enabled_approx_vertical(datasets, k=10, factor_k=2):

    partial = []
    for d in datasets:
        # d = np.float32(d)
        partial.append(local_SVD(d, k=factor_k * k).T)
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(partial, k=k, maxit=20, log=False)
    return H_i

def add_approx_vertical(datasets, k=10, factor_k=2):
    partial = []
    for d in datasets:
        # d = np.float32(d)
        partial.append(local_SVD(d.T, k=factor_k * k).T)
    H_i = np.sum(partial, axis=0)
    H_i, r = la.qr(H_i, mode='economic')
    g = [np.dot(d.T, H_i) for d in datasets]
    return g, H_i

def simulate_federated_horizontal_pca(datasets, k=10, factor_k=2):
    partial = []
    for d in datasets:
        #d = np.float32(d)
        partial.append(local_SVD(d, k=factor_k*k))
    print('Intermediate dimensions' +str(factor_k*k))
    dpca = aggregate_partial_SVDs(partial, factor_k*k)
    return dpca



def local_SVD(data, k=20):
    """
    Performs a singular value decomposition local data
    :param cov: A covariance matrix
    :param r: The number of top principal components to be considered
    :return: U_r*S_r (The product of the matrices taking the top r colums/rows)
    """

    # returns column vectors only
    U, S, UT, nd = sh.svd_sub(data, k)
    # In case we want to use more values in the approximation
    nd = min(nd, k)
    R = np.zeros((nd, nd))
    np.fill_diagonal(R, S[0:nd])
    U_r = UT[:, 0:nd]
    P = np.dot(np.sqrt(R), U_r.T)
    print(P.shape)
    return P

def aggregate_partial_SVDs(svds, t2=10):
    """
    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:
    """

    svds = np.concatenate(svds, axis=0)
    ndim = min(t2, svds.shape[0] - 1)

    print(svds.shape)
    U, S, UT, nd = sh.svd_sub(svds, ndim)
    return UT, S

def run_randomized(data_list, k, I,maxit, use_approximate=True, factor_k=2,filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None):
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(data_list,
                                                                              k=factor_k*k,
                                                                                maxit= I,
                                                                               filename=filename,
                                                                               u=u,
                                                                               choices=choices,
                                                                               precomputed_pca=precomputed_pca,
                                                                               fractev=fractev,
                                                                               federated_qr=federated_qr,
                                                                               v=v,
                                                                                gradient=gradient,
                                                                               epsilon=epsilon, log=log,
                                                                               g_ortho_freq=g_ortho_freq,
                                                                               g_init = g_init,
                                                                               previous_iterations=0)

    if use_approximate:
        for d in data_list:
            # d = np.float32(d)
            dl = local_SVD(d.T, k=factor_k * k).T
            H_stack.append(dl)
    H_stack = np.concatenate(H_stack, axis=1)
    H, S, G = lsa.svds(H_stack, k=H_stack.shape[1]-1)
    H = np.flip(H, axis=1)
    p = [np.dot(H.T,d) for d in data_list]
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list  = simulate_subspace_iteration(p,
                                                          k=k,
                                                            maxit= maxit,
                                                           filename=filename,
                                                           u=u,
                                                           choices=choices,
                                                           precomputed_pca=precomputed_pca,
                                                           fractev=fractev,
                                                           federated_qr=federated_qr,
                                                           v=None, # wrong vector
                                                            gradient=gradient,
                                                           epsilon=epsilon, log=log,
                                                           g_ortho_freq=g_ortho_freq,
                                                           g_init = G_i[:,0:k],
                                                           previous_iterations=iterations)
    G_i, R = la.qr(G_i, mode='economic')
    return G_i

def run_randomized_2(data_list, k, I,maxit, factor_k=2, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None):
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(data_list,
                                                                              k=k*factor_k,
                                                                                I= I,
                                                                               filename=filename,
                                                                               u=u,
                                                                               choices=choices,
                                                                               precomputed_pca=precomputed_pca,
                                                                               fractev=fractev,
                                                                               federated_qr=federated_qr,
                                                                               v=v,
                                                                                gradient=gradient,
                                                                               epsilon=epsilon, log=log,
                                                                               g_ortho_freq=g_ortho_freq,
                                                                               g_init = g_init,
                                                                               previous_iterations=1)

    H_stack = np.concatenate(H_stack, axis=1)
    H, S, G = lsa.svds(H_stack, k=H_stack.shape[1]-1)
    H = np.flip(H, axis=1)
    p = [np.dot(H.T,d) for d in data_list]
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(p,
                                                                           k=k,
                                                                           maxit=maxit,
                                                                           filename=filename,
                                                                           u=u,
                                                                           choices=choices,
                                                                           precomputed_pca=precomputed_pca,
                                                                           fractev=fractev,
                                                                           federated_qr=federated_qr,
                                                                           v=None,# wrong vector
                                                                           gradient=gradient,
                                                                           epsilon=epsilon, log=log,
                                                                           g_ortho_freq=g_ortho_freq,
                                                                           g_init=G_list,
                                                                           previous_iterations=iterations)
    return G_i





if __name__ == '__main__':
    if __name__ == '__main__':
        parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
        parser.add_argument('-f', metavar='file', type=str, help='filename of data file; default tab separated')
        parser.add_argument('--filetype', metavar='filetype', type=str, help='Type of the dataset')
        parser.add_argument('--sep', metavar='sep', type=str, help='spreadsheet separator, default tab', default='\t')
        parser.add_argument('--variance', action='store_true', help='spreadsheet separator, default tab')
        parser.add_argument('--center', action='store_true', help='center data')
        parser.add_argument('-o', metavar='outfile', type=str, help='output directory')
        parser.add_argument('-r', metavar='repeats', type=int, default=20, help='Number of times to repeat experiment')
        parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
        parser.add_argument('-t', metavar='tolerance', default=1e-9, type=float, help='Convergence tolerance')
        parser.add_argument('-s', metavar='sites', default='2,3,5,10', type=str,
                            help='comma separated list of number of sites to simulate, parsed as string')
        parser.add_argument('-i', metavar='iteration', default=2000, type=int, help='Maximum number of iterations')
        parser.add_argument('--header', metavar='iteration', default=None, type=int, help='header lines')
        parser.add_argument('--rownames', metavar='iteration', default=None, type=int, help='rownames')
        parser.add_argument('--names', metavar='iteration', default=None, type=str, help='names')
        parser.add_argument('--compare_pca', metavar='compare', default=None, type=str,
                            help='filename of precomputed pca to be compared to')
        parser.add_argument('--orthovector', metavar='compare', default=None, type=str,
                            help='filename of orthogonal file')
        parser.add_argument('--scaled', action='store_true', help='data is prescaled')
        parser.add_argument('--unequal', default=None, type=str, help='split unequal, load split file')
        parser.add_argument('--vert', action='store_true', help='run vertical split test')
        parser.add_argument('--hor', action='store_true', help='run horizontal split test')
        parser.add_argument('--ortho_freq', type=int, default=1, help='orthonormalisatio frequency for G')
        args = parser.parse_args()

        np.random.seed(95)
        # import scaled SNP file
        path = args.f
        filetype = args.filetype
        sep = args.sep
        k = args.k

        if args.names is None:
            dataset_name = os.path.basename(args.f)
        else:
            dataset_name = args.names

        if args.unequal is None:
            s = args.s
            splits = s.strip().split(',')
            splits = np.int8(splits)
            unequal = False
        else:
            unequal = True
            split_data = pd.read_csv(args.unequal, sep='\t', header=None)
            splits = []
            for i in range(split_data.shape[0]):
                l = split_data.iloc[i, :].tolist()
                cleanedList = [x for x in l if not np.isnan(x)]
                splits.append(cleanedList)

        maxit = args.i
        nr_repeats = args.r
        outdir = args.o
        scale = args.variance
        center = args.center
        ortho_freq = args.ortho_freq

        print(outdir)
        nr_samples = 0
        nr_features = 0
        if filetype == 'delim-list':
            data_list = []
            for f in path.split(','):
                data, sample_ids, variable_names = si.data_import(f, sep=sep)
                if scale or center:
                    data = si.scale_center_data_columnwise(data, center=center, scale_variance=scale)
                nr_samples += data.shape[0]
                nr_features += data.shape[1]
                data_list.append(data)
            data = data_list



        elif filetype == 'delim':
            data, sample_ids, variable_names = si.data_import(path, sep=sep, header=args.header, rownames=args.rownames)
            if scale or center:
                data = si.scale_center_data_columnwise(data, center=center, scale_variance=scale)
                nr_samples = data.shape[0]
                nr_features = data.shape[1]

        elif filetype == 'mnist':
            data, test_lables = mi.load_mnist(path, 'train')
            data = coo_matrix.asfptype(data)
            if scale or center:
                data = si.scale_center_data_columnwise(data, center=center, scale_variance=scale)
                nr_samples = data.shape[0]
                nr_features = data.shape[1]
            data = data.T

        elif filetype == 'gwas':
            bim = path + '.bim'
            traw = path + '.traw'

            if not args.scaled:
                traw_nosex = gi.remove_non_autosomes(bim, traw)

                data = gi.read_scale_write(infile=traw_nosex, outfile=path + '.traw.scaled', maf=0.01)

            else:
                data = pd.read_table(path + '.traw.scaled', header=None, sep='\t')
                data = data.values
            nr_samples = data.shape[0]
            nr_features = data.shape[1]
        else:
            raise Exception("Filetype not supported")

        if args.compare_pca is not None:
            precomputed_pca = pd.read_table(args.compare_pca, header=0, sep='\t')
            precomputed_pca = precomputed_pca.values
        else:
            precomputed_pca = None

        # vertical test
        if args.vert:
            vertical = op.join(outdir, 'vertical')
            os.makedirs(vertical, exist_ok=True)
            benchmark_vertical_approximate_pca(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                          outdir=vertical, precomputed_pca=precomputed_pca, unequal=unequal)

    # print('test')
    # start = time.monotonic()
    # import pnumpy as pn
    #
    # data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # # data, test_labels = mi.load_mnist(input_dir, 'train')
    # data = coo_matrix.asfptype(data)
    #
    # dataset_name = 'mnist'
    # maxit = 250
    # nr_repeats = 1
    # k = 10
    # splits = [5, 10]
    # outdir = '/home/anne/Documents/featurecloud/pca/approximative-vertical/results'
    # benchmark_vertical_approximate_pca(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-9,
    #                                        unequal=False, precomputed_pca=None)
    # print(time.monotonic()-start)