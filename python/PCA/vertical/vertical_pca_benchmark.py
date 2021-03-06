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
#import python.import_export.gwas_import as gi
import python.import_export.mnist_import as mi
import python.import_export.spreadsheet_import as si
import python.PCA.vertical.federated_qr as qr
import python.PCA.horizontal.power_iteration as powerit
import json
from python.PCA.horizontal.horizontal_pca_power_iteration import simulate_distributed_horizontal
from python.PCA.logging import *
from python.PCA.logging import TimerCounter

####### MATRIX POWER ITERATION SCHEME #######
def simulate_subspace_iteration(local_data, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, g_ortho_freq=1, g_init = None,
                                previous_iterations=None):
    """
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    """

    print('Orthonormalisation frequency'+ str(g_ortho_freq))
    G_list = []


    convergedH = False
    total_len = 0
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
        print(d.shape)
    start = 0

    if g_init is None:
        G_i = sh.generate_random_gaussian(total_len, k)
        G_i, R = la.qr(G_i, mode='economic')
        iterations = 0
    else:
        G_i = g_init
        iterations = 1

    if previous_iterations is not None:
        iterations = previous_iterations

    tol = TransmissionLogger()
    tol.open(filename)

    aol = AccuracyLogger()
    aol.open(filename)
    # send parts to local sites
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        tol.log_transmission( "G_i=SC", iterations, i, G_list[i])
        start = start + local_data[i].shape[1]

    # Initial guess
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)
    G_i_prev = G_i
    converged_eigenvals = []
    H_stack = []
    eigenvals_prev = None
    # Convergence can be reached when eigenvectors have converged, the maximal number of
    # iterations is reached or a predetermined number of eignevectors have converged.

    mot = TimerCounter()

    while not convergedH and iterations < maxit and len(converged_eigenvals) < k * fractev:
        iterations = iterations + 1
        #print(iterations)
        # add up the H matrices
        H_i = np.zeros((local_data[0].shape[0], k))
        for i in range(len(local_data)):
            # send local H matrices to server
            mot.start()
            H_local = np.dot(local_data[i], G_list[i])
            mot.stop()
            tol.log_transmission("H_local=CS", iterations, i, H_local)
            # add up H matrices at server and send them back to the clients
            H_i = H_i + H_local

        # Log only once for one site
        tol.log_transmission( "H_global=SC", iterations, 1, H_i)

        # free orthonormalisation in terms of transmission cost
        mot.start()
        H_i, R = la.qr(H_i, mode='economic')
        mot.stop()

        # Eigenvector update
        for i in range(len(G_list)):
            # Use gradient based update of the Eigenvectors

            if gradient:
                G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            else:
                # Use power iterations based update of the eigenvalue scheme
                mot.start()
                G_list[i] = np.dot(local_data[i].T, H_i)
                mot.stop()
                if not federated_qr:
                    tol.log_transmission("Gi_local=CS", iterations, i, G_list[i])

        # This is just for logging purposes
        G_i = np.concatenate(G_list, axis=0)

        # Eigenvalues are the norms of the eigenvecotrs
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))

        # this is not timed because it is done for logging purpose and not algorithmic reasons
        G_i, R = la.qr(G_i, mode='economic')


        convergedH, deltaH = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        # use guos convergence criterion for comparison
        #convergedH, deltaH = sh.convergence_checker_rayleigh(H_i, H_i_prev, eigenvals, eigenvals_prev, epsilon=1e-11)
        # just out of curiousity, log the
        #convergedG, deltaG = sh.eigenvector_convergence_checker(G_i, G_i_prev, tolerance=epsilon)
        H_i_prev = H_i
        G_i_prev = G_i
        if iterations < 10:
            H_stack.append(H_i)

        aol.log_current_accuracy(u=u, G_i=G_i, eigenvals=eigenvals, conv=deltaH, current_iteration=iterations,
                                 choices=choices, precomputed_pca=precomputed_pca, v=v, H_i=H_i)
    # log the time for matrix operations
#    log_time_keywords(filename, 'matrix_operations-subspace_iteration', mot.total())
    tol.close()
    aol.close()
    ortho, G_list = qr.simulate_federated_qr(G_list, encrypt=False)
    return G_i, eigenvals, converged_eigenvals, H_i, H_stack, iterations, G_list


####### ORIGINAL POWER ITERATION SCHEME #######
def residuals(V, a=None, sums=None):
    """
    Compute the residuals according to the formula
    """
    if a is None:
        a = sh.generate_random_gaussian(1, V.shape[0])
    sum = np.zeros(V.shape[0])
    if sums is not None:
        for v in range(V.shape[1]):
            sum = sum + sums[v] * V[:, v].T
    else:
        for v in range(V.shape[1]):
            sum = sum + np.dot(a, V[:, v:v + 1]) * V[:, v].T
    ap = a - sum
    return ap


def simulate_guo(local_data, maxit, V_k=None, starting_vector=None, filename=None, u=None, choices=None,
                 precomputed_pca=None, federated_qr=False, v=None, gradient=False, epsilon=1e-9, guo_epsilon=1e-11,
                 log=True, previous_iterations=None):
    """
    Retrieve the first eigenvector
    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of countermensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    """
    tol = TransmissionLogger()
    tol.open(filename)
    aol = AccuracyLogger()
    aol.open(filename)
    if previous_iterations is None:
        iterations = 0  # iteration counter
    else:
        iterations = previous_iterations

    # allow maxit for very eigenvector
    convergedH = False
    total_len = 0  # total number of samples/individuals in data
    for d in local_data:
        total_len = total_len + d.shape[1]
    if V_k is not None:
        id = V_k.shape[1]
    else:
        id = 0
    if starting_vector is None:
        # if there is no starting vector we generate an orthogonal
        # vector and start iterating
        if V_k is not None:
            # if it is not the first eigenvector use residuals
            # to orthogonalise
            G_i = residuals(V_k).T
        else:
            # if it is the first eigenvecto just generate
            # randonly and scale to unit norm
            G_i = sh.generate_random_gaussian(total_len, 1)
            G_i = G_i / np.linalg.norm(G_i)  # normalize
    else:
        # the vector is already preiterated and is orthonormal, we
        # just have to assure it is a 2d array.
        G_i = np.reshape(starting_vector, (total_len, 1))

    start = 0  # running variable to partition G_i
    G_list = []  # this are the parital eigevenctors
    Vk_list = []  # these are the partial eigenvectors already iterated
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        tol.log_transmission("G_i=SC", iterations, i, G_list[i], id+1)
        if V_k is not None:
            Vk_list.append(V_k[start:start + local_data[i].shape[1], :])
        start = start + local_data[i].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)
    G_i_prev = G_i
    gi_norm_prev = la.norm(H_i_prev)
    mot = TimerCounter()
    while not convergedH and iterations < maxit:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], 1))  # dummy initialise the H_i matrix
        mot.start()
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            tol.log_transmission( "H_local=CS", iterations, i, H_local ,id+1)
            H_i = H_i + H_local
        tol.log_transmission( "H_global=SC", iterations, 1, H_i ,id+1)

        for i in range(len(local_data)):
            if gradient:
                G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            else:
                G_list[i] = np.dot(local_data[i].T, H_i)
        mot.stop()
        gi_norm = 0
        if V_k is None:
            # compute the norm of the eigenvector and done.
            for i in range(len(G_list)):
                local_norm = np.sum(np.square(G_list[i]))
                gi_norm = gi_norm + np.sum(np.square(G_list[i]))
                tol.log_transmission( "local_norm=CS", iterations, i, local_norm ,id+1)

        elif federated_qr:
            temp = []
            for i in range(len(G_list)):
                gi_norm = gi_norm + np.sum(np.square(G_list[i]))
                temp.append(np.concatenate([Vk_list[i], G_list[i]], axis=1))

            G_i, G_list = qr.simulate_federated_qr(temp, encrypt=False)

            for i in range(len(G_list)):
                G_list[i] = G_list[i][:, id:id + 1]

        else:
            local_sums = []
            for i in range(len(G_list)):
                sum = []
                for vi in range(Vk_list[i].shape[1]):
                    dp = np.dot(G_list[i].T, Vk_list[i][:, vi:vi + 1]).flatten()
                    sum.append(dp)
                # cast to numpy array to determine size
                tol.log_transmission( "local_dot_prod=CS", iterations, i, np.asarray(sum), id+1)
                local_sums.append(sum)
            local_sums = np.asarray(local_sums)
            local_sums = np.sum(local_sums, axis=0).flatten()  # flatten to remove nesting
            tol.log_transmission( "global_dot_prod=SC", iterations, 1, local_sums, id+1)

            for i in range(len(G_list)):
                ap = G_list[i]
                for vi in range(Vk_list[i].shape[1]):
                    it = local_sums[vi] * Vk_list[i][:, vi:vi + 1].T
                    it = np.reshape(it, ap.shape)
                    ap = ap - it
                G_list[i] = ap

            for i in range(len(G_list)):
                c_n = np.sum(np.square(G_list[i]))
                tol.log_transmission("local_norm=CS", iterations, i, c_n, id+1)
                gi_norm = gi_norm + c_n

        gi_norm = np.sqrt(gi_norm)
        if V_k is None or not federated_qr:
            for i in range(len(G_list)):
                G_list[i] = G_list[i] / gi_norm
                if log:
                    # subsequent orthonormalisation at agrgegator
                    # (only current G_i because rest is assumed to be stored at agrgegator)
                    tol.log_transmission( "ALT_G_i_local=CS", iterations, i, G_list[i], id+1)
                    tol.log_transmission( "ALT_G_i=SC", iterations, i, G_list[i], id+1)


        #if gradient:
        convergedH, deltaH = sh.convergence_checker_rayleigh(H_i, H_i_prev, [gi_norm], [gi_norm_prev] ,epsilon=guo_epsilon)
        #else:
        #    convergedH, deltaH = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        #    print(convergedH)
        #    print(deltaH)
        H_i_prev = H_i
        gi_norm_prev = gi_norm
        tol.log_transmission( "global_norm=SC", iterations, 1, gi_norm, id+1)

        G_i = np.concatenate(G_list, axis=0)
        convergedG, deltaG = sh.eigenvector_convergence_checker(G_i, G_i_prev, tolerance=epsilon)
        G_i_prev = G_i
        aol.log_current_accuracy(u[:, id:id + 1], G_i, eigenvals=[gi_norm], conv=deltaH, current_iteration=iterations,
                             choices=choices, precomputed_pca=precomputed_pca, current_ev=id + 1,
                                v=v[:, id:id + 1], H_i=H_i)
    # return a complete eigenvector
    #print(gi_norm)
    print(iterations)
    log_time_keywords(filename, 'matrix_operations-subspace_iteration', mot.total())
    tol.close()
    aol.close()
    return G_i, iterations


def compute_k_eigenvectors(data_list, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None,
                           federated_qr=False, v=None, gradient=True, epsilon=1e-9, guo_epsilon=1e-11):
    ug, it = simulate_guo(data_list, maxit=maxit, filename=filename, u=u, choices=choices, federated_qr=federated_qr, v=v,
                      gradient=gradient, epsilon=epsilon, guo_epsilon=guo_epsilon)
    u_all = ug
    for i in range(1, k):
        ug2, it= simulate_guo(data_list, maxit=maxit+it, V_k=u_all, filename=filename, u=u, choices=choices,
                           federated_qr=federated_qr, v=v, gradient=gradient, epsilon=epsilon, guo_epsilon=guo_epsilon,
                              previous_iterations=it)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all




####### BENCHMARK RUNNER #######
def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-9, precomputed_pca=None,
                  unequal=False, horizontal=False, guo_epsilon=1e-11, ortho_freq=1):
    """
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

    """
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

            # run power iteration benchmark
            # simultaneous fully federated QR
            grad = False
            grad_name = 'power'
            mode = 'federated_qr'
            fedqr = True
            print('power - matrix - '+ mode)
            outdir_gradient = op.join(outdir, 'matrix', str(s), grad_name, mode, str(1))
            os.makedirs(outdir_gradient, exist_ok=True)
            filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                   precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=1)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)

            # simultaneous only H
            grad = False
            grad_name = 'power'
            print('power - matrix - ' + mode)
            outdir_gradient = op.join(outdir, 'matrix', str(s), grad_name, mode, str(ortho_freq))
            os.makedirs(outdir_gradient, exist_ok=True)
            filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
            simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                        precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                        epsilon=epsilon, g_ortho_freq=ortho_freq)
            end = time.monotonic()
            log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)

            # Run power iteration based benchmark
            # Sequential
            #print('power - sequential - '+ mode)
            #outdir_gradient = op.join(outdir, 'vector', str(s), grad_name, mode)
            #os.makedirs(outdir_gradient, exist_ok=True)
            #filename = create_filename(outdir_gradient, dataset_name_guo + '_' + mode, s, c, k, maxit, start)

            #start = time.monotonic()
            #compute_k_eigenvectors(data_list, k=k, maxit=maxit, u=u, filename=filename, choices=choice,
            #                      precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
            #                       epsilon=epsilon, guo_epsilon=guo_epsilon)
            #end = time.monotonic()
            #log_time(logftime, 'guo_single' + '_' + mode, end - start, s, c)

            # Run Guo version
            # Sequention
            grad = True
            fedqr = False
            grad_name = 'gradient'
            mode = 'central_qr'
            print('gradient - sequential - '+ mode)
            outdir_gradient = op.join(outdir, 'vector', str(s), grad_name, mode, str(1))
            os.makedirs(outdir_gradient, exist_ok=True)

            filename = create_filename(outdir_gradient, dataset_name_guo + '_' + mode, s, c, k, maxit, start)

            start = time.monotonic()
            compute_k_eigenvectors(data_list, k=k, maxit=maxit, u=u, filename=filename, choices=choice,
                                   precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                   epsilon=epsilon, guo_epsilon=guo_epsilon)
            end = time.monotonic()
            log_time(logftime, 'guo_single' + '_' + mode, end - start, s, c)

            logf = op.join(outdir, 'log_choices.log')
            log_choices(logf, filename, choice)

            if horizontal:
                start = time.monotonic()
                outdir_pow = op.join(outdir, 'power_iteration', str(s))
                os.makedirs(outdir_pow, exist_ok=True)
                filename = create_filename(outdir_pow, dataset_name, s, c, k, maxit, start)
                simulate_distributed_horizontal(data_list, p=k, v=v, filename=filename, choices=choice, maxit=maxit)
                end = time.monotonic()
                log_time(logftime, 'powerit_hor', end - start, s, c)


####### BENCHMARK RUNNER #######

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
    parser.add_argument('--ortho_freq',type=int, default=1, help='orthonormalisatio frequency for G')
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
    ortho_freq=args.ortho_freq

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
        the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                      outdir=vertical, precomputed_pca=precomputed_pca, unequal=unequal, ortho_freq=ortho_freq)

    # horizontal test
    if args.hor:
        horizontal = op.join(outdir, 'horizontal')
        os.makedirs(horizontal, exist_ok=True)
        the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                      outdir=horizontal, precomputed_pca=precomputed_pca, unequal=unequal, horizontal=True)
