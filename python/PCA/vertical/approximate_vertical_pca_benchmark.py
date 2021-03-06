from python.PCA.vertical.vertical_pca_benchmark import simulate_subspace_iteration, simulate_guo, compute_k_eigenvectors
from python.evaluation.data_aggregation import create_dataframe
import python.PCA.logging as l
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


import python.PCA.shared_functions as sh
import python.import_export.mnist_import as mi
import python.import_export.spreadsheet_import as si
#import python.import_export.gwas_import as gi

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
    dataset_name_guo = dataset_name + '_guo'
    guo_epsilon = 1e-11
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
            ortho_freq = maxit+1 # will not be reached
            fedqr = False

            # simulate the run
            start = time.monotonic()
            mode = 'randomized-no-approx'
            outdir_approx = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            os.makedirs(outdir_approx, exist_ok=True)
            filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)

            run_randomized(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
                           precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                           epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None, use_approximate=False)
            end = time.monotonic()
            log_time(logftime, mode, end - start, s, c)
            print(mode + ' ' + str(end - start))

            # start = time.monotonic()
            # mode = 'randomized-approx-projected'
            # outdir_approx = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            # os.makedirs(outdir_approx, exist_ok=True)
            # filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            #
            # run_randomized_2(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
            #                precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None)
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))


            #
            #
            # start = time.monotonic()
            # mode = 'randomized-approx-initial'
            # outdir_approx = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            # os.makedirs(outdir_approx, exist_ok=True)
            # filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            #
            # run_randomized(data_list, k, I=10, maxit=maxit, u=u, filename=filename, choices=choice,
            #                precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None, use_approximate=True)
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))
            #
            #
            # start = time.monotonic()
            # mode = 'random-init'
            # print('power - matrix - ' + mode)
            # outdir_gradient = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            # os.makedirs(outdir_gradient, exist_ok=True)
            # filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
            # simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
            #                             precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                             epsilon=epsilon, g_ortho_freq=ortho_freq)
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))
            #
            #
            # # Compute federated approximate PCA
            # start = time.monotonic()
            # mode = 'approximative-init'
            # outdir_approx = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            # os.makedirs(outdir_approx, exist_ok=True)
            # filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            # mot =TimerCounter()
            # mot.start()
            # aol = AccuracyLogger()
            # aol.open(filename)
            # g, h = approximate_vertical(data_list, k, factor_k=2)
            # mot.stop()
            # log_time_keywords(filename, 'matrix_operations-approximative-init', mot.total())
            # g = np.concatenate(g, axis=0)
            # g,r  = la.qr(g, mode='economic')
            # aol.log_current_accuracy(u=u, G_i=g, current_iteration=1, precomputed_pca=precomputed_pca, v=v, H_i=h, choices= choice)
            # end = time.monotonic()
            # aol.close()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))
            #
            # #Compute full decomposition using approximative PCA as seed
            # start = time.monotonic()
            # simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
            #                             precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                             epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k])
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))
            #
            #
            # # # Compute federated approximate PCA
            # start = time.monotonic()
            # mode = 'approximative-init-smpc'
            # outdir_approx = op.join(outdir, 'matrix', str(s), mode, str(ortho_freq))
            # os.makedirs(outdir_approx, exist_ok=True)
            # filename = create_filename(outdir_approx, dataset_name + '_' + mode, s, c, k, maxit, start)
            # mot =TimerCounter()
            # mot.start()
            # aol.open(filename)
            # g, h = approximate_vertical_smpc(data_list, k, factor_k=2, u=u, filename=filename, choices=choice,
            #                             precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                             epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=None)
            # mot.stop()
            # log_time_keywords(filename, 'matrix_operations-approximative-init-smpc', mot.total())
            # g = np.concatenate(g, axis=0)
            # aol.log_current_accuracy(u=u, G_i=g, current_iteration=1,
            #                       precomputed_pca=precomputed_pca, v=v, H_i=h, choices=choice)
            # end = time.monotonic()
            # aol.close()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode+' '+str(end-start))
            #
            # # Compute full decomposition using approximative PCA as seed
            # start = time.monotonic()
            # # iterations are hardcoded.
            # simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
            #                             precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                             epsilon=epsilon, g_ortho_freq=ortho_freq, g_init=g[:, 0:k], previous_iterations=20)
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)
            # print(mode + ' ' + str(end - start))
            #
            # # Run Guo version
            # # Sequention
            # grad = True
            # fedqr = False
            # mode = 'guo'
            # print('gradient - sequential - ' + mode)
            # outdir_gradient = op.join(outdir, 'vector', str(s), mode, str(1))
            # os.makedirs(outdir_gradient, exist_ok=True)
            # filename = create_filename(outdir_gradient, dataset_name_guo + '_' + mode, s, c, k, maxit, start)
            #
            # start = time.monotonic()
            # compute_k_eigenvectors(data_list, k=k, maxit=maxit, u=u, filename=filename, choices=choice,
            #                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                        epsilon=epsilon, guo_epsilon=guo_epsilon)
            # end = time.monotonic()
            # log_time(logftime, mode, end - start, s, c)


def approximate_vertical(data_list, k=10, factor_k=2, filename=None):
    data_list = [d.T for d in data_list]
    v, e = simulate_federated_horizontal_pca(data_list, k=k, factor_k=factor_k, filename=filename)
    g = [np.dot(d, v) for d in data_list]
    return g, v

def approximate_vertical_smpc(data_list, k=10, factor_k=2,filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None):
    data_list = [d.T for d in data_list]
    v = smpc_enabled_approx_vertical(data_list, k=k, factor_k=factor_k,filename=filename,
                                                u=u,
                                                choices=choices,
                                                precomputed_pca=precomputed_pca,
                                                fractev=fractev,
                                                federated_qr=federated_qr,
                                                v=v,
                                                gradient=gradient,
                                                epsilon=epsilon,
                                                log=log,
                                                g_ortho_freq=g_ortho_freq,
                                                g_init=g_init)
    g = [np.dot(d, v) for d in data_list]
    return g, v

def smpc_enabled_approx_vertical(datasets, k=10, factor_k=2,filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None):

    partial = []
    for d in datasets:
        # d = np.float32(d)
        partial.append(local_SVD(d, k=factor_k * k).T)
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(partial, k=k, maxit=20,
                                                                                                        filename=filename,
                                                                                                        u=None,
                                                                                                        choices=choices,
                                                                                                        precomputed_pca=precomputed_pca,
                                                                                                        fractev=fractev,
                                                                                                        federated_qr=federated_qr,
                                                                                                        v=None,
                                                                                                        gradient=gradient,
                                                                                                        epsilon=epsilon,
                                                                                                        log=log,
                                                                                                        g_ortho_freq=g_ortho_freq,
                                                                                                        g_init=g_init,
                                                                                            previous_iterations=0)
    return H_i

def simulate_federated_horizontal_pca(datasets, k=10, factor_k=2, filename=None):
    partial = []
    tol = TransmissionLogger()
    tol.open(filename)
    i = 0
    for d in datasets:
        #d = np.float32(d)
        dl = local_SVD(d, k=factor_k*k)
        partial.append(dl)
        tol.log_transmission( "H_local=CS", -1, i, dl)
        i =i+1
    print('Intermediate dimensions' +str(factor_k*k))
    dpca = aggregate_partial_SVDs(partial, factor_k*k)
    tol.log_transmission( "H_global=SC", -2, 1, dpca)
    tol.close()
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
    mot =TimerCounter()
    tol = TransmissionLogger()
    tol.open(filename)
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
    mot.start()

    if use_approximate:
        i = 0
        for d in data_list:
            # d = np.float32(d)
            dl = local_SVD(d.T, k=factor_k * k).T
            tol.log_transmission( "H_local=CS", iterations+1, i, dl)
            H_stack.append(dl)
            i = i+1
    H_stack = np.asarray(np.concatenate(H_stack, axis=1))
    H, S, G = lsa.svds(H_stack, k=H_stack.shape[1]-1)
    tol.log_transmission( "H_global=SC", iterations+1, 1, H)
    mot.stop()
    tol.close()
    H = np.flip(H, axis=1)
    p = [np.dot(H.T,d) for d in data_list]
    covs = [np.dot(p1, p1.T) for p1 in p]
    u1,s1, v1 = lsa.svds(np.sum(covs, axis=0), k=k)
    u1 = np.flip(u1, axis=1)
    g1 = [np.dot(p1.T, u1) for p1 in p]
    G_i = np.concatenate(g1, axis=0)
    # G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list  = simulate_subspace_iteration(p,
    #                                                       k=k,
    #                                                         maxit= maxit,
    #                                                        filename=filename,
    #                                                        u=u,
    #                                                        choices=choices,
    #                                                        precomputed_pca=precomputed_pca,
    #                                                        fractev=fractev,
    #                                                        federated_qr=federated_qr,
    #                                                        v=None, # wrong vector
    #                                                         gradient=gradient,
    #                                                        epsilon=epsilon, log=log,
    #                                                        g_ortho_freq=g_ortho_freq,
    #                                                        g_init = G_i[:,0:k],
    #                                                        previous_iterations=iterations)

    G_i, R = la.qr(G_i, mode='economic')
    aol = AccuracyLogger()
    aol.open(filename)
    aol.log_current_accuracy(u=u, G_i=G_i, eigenvals=eigenvals, conv=None, current_iteration=iterations,
                             choices=choices, precomputed_pca=precomputed_pca, v=v, H_i=H_i)
    aol.close()
    #log_time_keywords(filename, 'matrix_operations-randomized', mot.total())
    return G_i

def run_randomized_2(data_list, k, I,maxit, factor_k=2, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None,
                            use_approximate=True):
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(data_list,
                                                                              k=k*factor_k,
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
    mot = TimerCounter()
    mot.start()
    H_stack = np.concatenate(H_stack, axis=1)
    H, S, G = lsa.svds(H_stack, k=H_stack.shape[1]-1)
    H = np.flip(H, axis=1)
    p = [np.dot(H.T,d) for d in data_list]
    if use_approximate:
        G_i,vv = approximate_vertical(p, k=k, factor_k=2, filename=filename)
        G_i = np.concatenate(G_i, axis=0)
    mot.stop()
    log_time_keywords(filename, 'matrix_operations-randomized', mot.total())
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
                                                                           g_init=G_i[:, 0:k],
                                                                           previous_iterations=iterations)
    #G_i, R = la.qr(G_i, mode='economic')
    G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list = simulate_subspace_iteration(data_list,
                                                                              k=k,
                                                                                maxit= 1+iterations,
                                                                               filename=filename,
                                                                               u=u,
                                                                               choices=choices,
                                                                               precomputed_pca=precomputed_pca,
                                                                               fractev=fractev,
                                                                               federated_qr=federated_qr,
                                                                               v=None,
                                                                                gradient=gradient,
                                                                               epsilon=epsilon, log=log,
                                                                               g_ortho_freq=g_ortho_freq,
                                                                               g_init = G_i[:, 0:k],
                                                                               previous_iterations=iterations)
    print(co.compute_angles(H_i, v))
    return G_i





if __name__ == '__main__':
    local=True
    if local:
        start = time.monotonic()

        data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
        # data, test_labels = mi.load_mnist(input_dir, 'train')
        data = coo_matrix.asfptype(data)

        dataset_name = 'mnist'
        maxit = 1000
        nr_repeats = 1
        k = 10
        splits = [5, 10]
        outdir = '/home/anne/Documents/featurecloud/pca/approximative-vertical/results1'
        benchmark_vertical_approximate_pca(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-9,
                                           unequal=False, precomputed_pca=None)
        print(time.monotonic() - start)
        outd = ['matrix', 'vector']
        for od in outd:
            basepath = op.join(outdir,od)
            create_dataframe(basepath=basepath, suffix='.angles.u')
            create_dataframe(basepath=basepath, suffix='.angles.v')
            create_dataframe(basepath=basepath, suffix='.mev.u')
            create_dataframe(basepath=basepath, suffix='.mev.v')
            create_dataframe(basepath=basepath, suffix='.transmission')
            create_dataframe(basepath=basepath, suffix='.eo')


    else:
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
        parser.add_argument('--sf', default=None, type=str, help='suffix scaled file')
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

                data = gi.read_scale_write_pandas(infile=traw_nosex, outfile=path + '.traw.scaled', maf=0.01, major_2=False)

            else:
                print('Loading scaled data')
                if args.sf is None:
                    data = pd.read_table(path + '.traw.scaled', header=None, sep='\t')
                else:
                    data = pd.read_table(path+ '.' +args.sf, header=None, sep='\t')
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

            outd = ['matrix', 'vector']
            for od in outd:
                basepath = op.join(vertical,od)
                create_dataframe(basepath=basepath, suffix='.angles.u')
                create_dataframe(basepath=basepath, suffix='.angles.v')
                create_dataframe(basepath=basepath, suffix='.mev.u')
                create_dataframe(basepath=basepath, suffix='.mev.v')
                create_dataframe(basepath=basepath, suffix='.transmission')
                create_dataframe(basepath=basepath, suffix='.eo')

