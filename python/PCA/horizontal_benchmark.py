import python.PCA.power_iteration_runner as por
import python.PCA.shared_functions as sh
import python.import_export.mnist_import as mi
from scipy.sparse import coo_matrix
import numpy as np
import scipy.sparse.linalg as lsa
import python.PCA.comparison as co
import argparse as ap
import os.path as path
import time
import numpy as np
import scipy.sparse.linalg as lsa
import scipy.linalg as la
import os.path as op

#import python.PCA.import_export.
import python.PCA.power_iteration as powerit
import python.PCA.proxy_covariance as dpca
import python.PCA.vertical_pca_library as vert
import python.PCA.comparison as co
import python.PCA.vertical_pca_benchmark as benchmark
import python.PCA.federated_qr as qr
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
import python.PCA.federated_qr as qr

def simulate_distributed(local, p=10, weights=None, tolerance=0.000001, scipy = None, maxit=2000, filename=None, choices=None, log = True):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[1]
    X = powerit.generate_random_gaussian(d, p, sigma=1)
    X, R = la.qr(X, mode='economic')
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    X_prev = X
    while (not converged) and count<maxit:
        count = count + 1
        locals = []
        for l in local:
            d1 = np.dot(l, X)
            d1 = np.dot(l.T, d1)
            locals.append(d1)
        X, E, converged = powerit.pooling_step(locals, X, weights=weights)
        converged3, conv, converged_eigenvals, delta  = vert.convergence_checker(X, X_prev, return_converged=True, epsilon=0.000001)
        if log:
            print(scipy.shape)
            print(X.shape)
            benchmark.log_current_accuracy(scipy=scipy, G_i=X, eigenvals=E, conv=delta, current_iteration=count,
                                 filename=filename, choices=range(scipy.shape[0]))
    ord = np.argsort(E)
    X = np.flip(X[:, ord], axis=1)
    E = np.flip(np.sort(E))
    print(local[0].shape)
    print(X.shape)
    return X, E, count

def simulate_guo_benchmark(local_data, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0, federated_qr = False, v = None):
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
    G_i_prev = G_i
    converged_eigenvals = []
    eigenvals_prev=None
    convergeda = False
    while not converged and iterations < maxit and len(converged_eigenvals)<k*fractev:
        iterations = iterations + 1
        #print(iterations)
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
        if not eigenvals_prev is None:
            convergeda, conva, converged_eigenvalsa, deltaa = gv.convergence_checker_a(H_i, H_i_prev, return_converged=True,alpha_norm=eigenvals, alpha_norm_prev=eigenvals_prev, epsilon = 0.000000001)
            #print('a'+ str(convergeda))
            #print(conva)
            #print(deltaa)
        eigenvals_prev = eigenvals
        eigenvals = np.sqrt(eigenvals)

        if not federated_qr:
            G_i, R = la.qr(G_i, mode='economic')
            start = 0
            for i in range(len(G_list)):
                G_list[i] = G_i[start:start + local_data[i].shape[1], :]
                #log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
                start = start + local_data[i].shape[1]
        else:
            G_i, G_list = qr.simulate_federated_qr(G_list, encrypt=False, filename=filename, repeat= iterations )


        converged, conv, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True, epsilon=0.01)
        #print('converged'+str(converged))
        gi_delta_obj = sh.eigenvector_convergence_checker(G_i, G_i_prev)
        H_i_prev = H_i
        G_i_prev  = G_i
        benchmark.log_current_accuracy(scipy=v, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations, filename=filename, choices=choices, precomputed_pca=precomputed_pca, gi_delta_obj = gi_delta_obj, v = u, H_i=H_i)
    G_i = np.concatenate(G_list)
    print(iterations)
    return G_i, eigenvals, converged_eigenvals, H_i

def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir, convergence_eps=1e-6, precomputed_pca=None, unequal=False, nr_samples=0, nr_features=0):
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

    current_split = 0
    for c in range(nr_repeats):
        for s in splits:
            # filename will be the same for angle log file and correlation log file
            timer = time.monotonic()


            # split the data
            if not islist:
                if unequal:
                    # if unequal then the number of sites is the length of s and s itself contains the splits
                    data_list, choice = sh.partition_data_horizontally(data, len(s), randomize=True, perc=s, equal=False)
                    s = current_split
                    current_split += 1
                else:
                    data_list, choice = sh.partition_data_horizontally(data, s, randomize=True)
            else:
                choice = range(data.shape[1])




            logftime = op.join(outdir, 'time.log')

            # # simulate the run

            start = time.monotonic()
            filename = benchmark.init_benchmark(outdir=outdir, dataset_name=dataset_name, maxit=maxit, counter=c, nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps, splits=s, timer=timer, transmission_costs=True)

            #simulate_distributed(data_list, p = k, maxit=maxit, scipy=u, filename=filename, choices=choice)
            end = time.monotonic()
            benchmark.log_time(logftime, 'naive_distributed', end - start, s, c)

            for d in range(len(data_list)):
                data_list[d] = data_list[d].T
            start = time.monotonic()

            filename2 = benchmark.init_benchmark(outdir=outdir, dataset_name=dataset_name+'_vert', maxit=maxit, counter=c,
                                                nr_samples=nr_samples, nr_features=nr_features, k=k,
                                                convergence_eps=convergence_eps, splits=s, timer=timer,
                                                transmission_costs=True)
            simulate_guo_benchmark(data_list, k, maxit=maxit, u=u, filename=filename2, choices=choice,
                                   precomputed_pca=precomputed_pca, federated_qr=False, v=v)
            end = time.monotonic()
            benchmark.log_time(logftime, 'smart_distributed', end - start, s, c)


if __name__ == '__main__':
    # data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # #
    # # data, sample_ids, variable_names = si.data_import('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t', header=0, rownames=0)
    # # data = si.scale_center_data_columnwise(data)
    # data = coo_matrix.asfptype(data)
    # outdir ='/home/anne/Documents/featurecloud/pca/vertical-pca/results_hori/mnist/unequal'
    # the_epic_loop(data, 'mnist', maxit=200, nr_repeats=20, k=10, splits=[2,3,5,10], outdir=outdir)
    #
    # u, s, v = lsa.svds(data, k=10)
    # # u = np.flip(u, axis=1)
    # # s = np.flip(s)
    # # v = np.flip(v.T, axis=1)
    # #
    # # datasplit, choice = sh.partition_data_horizontally(data=data, splits=3)
    # # x, e, count = por.simulate_distributed(datasplit, p=10, scipy=v)
    # #
    # # print(co.compute_angles(x, v))
    # #

    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; default tab separated')
    parser.add_argument('--filetype', metavar='filetype', type=str, help='Type of the dataset')
    parser.add_argument('--sep', metavar='sep', type=str, help='spreadsheet separator, default tab', default='\t')
    parser.add_argument('--variance', action='store_true', help='spreadsheet separator, default tab')
    parser.add_argument('--center', action='store_true', help='center data')
    parser.add_argument('-o', metavar='outfile', type=str, help='output directory')
    parser.add_argument('-r', metavar='repeats', type=int, default=20, help='Number of times to repeat experiment')
    parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
    parser.add_argument('-t', metavar='tolerance', default=1e-6, type=float, help='Convergence tolerance')
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
    args = parser.parse_args()

    np.random.seed(95)
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
            l =split_data.iloc[i, :].tolist()
            cleanedList = [x for x in l if not np.isnan(x)]
            splits.append(cleanedList)

    maxit = args.i
    nr_repeats = args.r
    outdir = args.o
    scale = args.variance
    center = args.center

    print(outdir)
    nr_samples = 0
    nr_features = 0

    if filetype == 'delim':
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

    the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                  outdir=outdir, precomputed_pca=precomputed_pca, unequal=unequal)
