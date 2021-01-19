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
import python.import_export.gwas_import as gi
import python.import_export.mnist_import as mi
import python.import_export.spreadsheet_import as si
import python.PCA.vertical.federated_qr as qr
import python.PCA.horizontal.power_iteration as powerit


####### LOGGING FUNCTIONS #######
def create_filename(outdir, dataset_name, splits, counter, k, maxit, time):
    """
    Make a file name
    Args:
        dataset_name: Name of the dataset currently running
        splits: The current split
        counter: Current repeat experiment
        k: Number of targeted dimensions
        maxit: maximal iterations

    Returns: A concatenated filename with filestamp

    """
    fn = dataset_name + '_' + str(splits) + '_' + str(counter) + '_' + str(k) + '_' + str(maxit) + '_' + str(time)
    fn = op.join(outdir, fn)
    return fn


def log_current_accuracy(u, G_i, eigenvals, conv, current_iteration, filename, choices, precomputed_pca=None,
                         current_ev=None, gi_delta_obj=None, v=None, H_i=None):
    """
    Log the current iterations angle to the canonical
    Args:
        u: column vector matrix with canonical eigenvectors
        G_i: column vector based matrix with current eigenvector estimation
        current_iteration: iteration index
        filename: output filename prefix > out will be saved to out.angles, and out.cor

    Returns: None

    """
    if current_ev is not None:
        info_string = str(current_ev) + '\t' + str(current_iteration)
    else:
        info_string = str(current_iteration)

    if not u is None:
        with open(filename + '.angles.u', 'a+') as handle:
            angles = co.compute_angles(u[choices, :], G_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

        with open(filename + '.cor', 'a+') as handle:
            correlations = co.compute_correlations(u[choices, :], G_i)
            if correlations is not None and len(correlations) > 0:
                info = cv.collapse_array_to_string(correlations, info_string)
                handle.write(info)

    if not v is None:
        with open(filename + '.angles.v', 'a+') as handle:
            angles = co.compute_angles(v, H_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

    with open(filename + '.eigenval', 'a+') as handle:
        info = cv.collapse_array_to_string(eigenvals, info_string)
        handle.write(info)

    if conv is not None:
        with open(filename + '.conv', 'a+') as handle:
            conv = cv.collapse_array_to_string(conv, info_string)
            handle.write(conv)

    if precomputed_pca is not None:
        with open(filename + '.angles_precomp', 'a+') as handle:
            angles = co.compute_angles(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(angles, info_string)
            handle.write(info)

        with open(filename + '.cor_precomp', 'a+') as handle:
            correlations = co.compute_correlations(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(correlations, info_string)
            handle.write(info)
    if gi_delta_obj is not None:
        with open(filename + '.eigenvector_convergence', 'a+') as handle:
            conv = cv.collapse_array_to_string(gi_delta_obj, info_string)
            handle.write(conv)

def log_choices(logfile, filename, choices):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(cv.collapse_array_to_string(choices, filename))


def log_transmission(logfile, log_entry_string, iterations, counter, element, eigenvector=10):
    with open(logfile + '.transmission', 'a+') as handle:
        if type(element) == 'numpy.ndarray':
            try:
                handle.write(log_entry_string + '\t' + str(iterations)
                             + '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(
                    sys.getsizeof(element.tobytes())) + '\n')
            except AttributeError:
                handle.write(log_entry_string + '\t' + str(iterations)
                             + '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(
                    sys.getsizeof(element)) + '\n')
        else:
            handle.write(log_entry_string + '\t' + str(iterations)
                         + '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(
                sys.getsizeof(element)) + '\n')


def log_time(logfile, algorithm, time, split, repeat):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(algorithm + '\t' + str(split) + '\t' + str(repeat) + '\t' + str(time) + '\n')


####### MATRIX POWER ITERATION SCHEME #######
def simulate_subspace_iteration(local_data, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=0.000001):
    """
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    """

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
        # log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
        start = start + local_data[i].shape[1]

    # Initial guess
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)
    G_i_prev = G_i
    converged_eigenvals = []
    eigenvals_prev = None
    # Convergence can be reached when eigenvectors have converged, the maximal number of
    # iterations is reached or a predetermined number of eignevectors have converged.
    while not converged and iterations < maxit and len(converged_eigenvals) < k * fractev:
        iterations = iterations + 1

        # add up the H matrices
        H_i = np.zeros((local_data[0].shape[0], k))
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            # log_transmission(filename, "H_local=CS", iterations, i, H_local)
            H_i = H_i + H_local
        # log_transmission(filename, "H_global=SC", iterations, 1, H_i)

        # free orthonormalisation in terms of transmission cost
        H_i, R = la.qr(H_i, mode='economic')

        # Eigenvector update
        for i in range(len(G_list)):
            # Use gradient based update of the Eigenvectors
            if gradient:
                G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            else:
                # Use power iterations based update of the eigenvalue scheme
                G_list[i] = np.dot(local_data[i].T, H_i)
            # log_transmission(filename, "Gi_local=CS", iterations, i, G_list[i])

        # This is just for logging purposes
        G_i = np.concatenate(G_list, axis=0)

        # Eigenvalues are the norms of the eigenvecotrs
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        # Save eigenvalues
        eigenvals = np.sqrt(eigenvals)
        eigenvals_prev = eigenvals

        if not federated_qr:
            # Centralised QR, just send the eigenvectors to the aggregator and orthonormalise
            # Concatenation happend
            G_i, R = la.qr(G_i, mode='economic')
            start = 0
            # redistribute
            for i in range(len(G_list)):
                G_list[i] = G_i[start:start + local_data[i].shape[1], :]
                # log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
                start = start + local_data[i].shape[1]
        else:
            G_i, G_list = qr.simulate_federated_qr(G_list, encrypt=False, filename=filename, repeat=iterations)

        convergedH, deltaH = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        # just out of curiousity, log the
        convergedG, deltaG = sh.eigenvector_convergence_checker(G_i, G_i_prev, tolerance=epsilon)
        H_i_prev = H_i
        G_i_prev = G_i
        log_current_accuracy(u=u, G_i=G_i, eigenvals=eigenvals, conv=deltaH, current_iteration=iterations,
                             filename=filename, choices=choices, precomputed_pca=precomputed_pca,
                             gi_delta_obj=deltaG, v=v, H_i=H_i)
    G_i = np.concatenate(G_list)
    print(iterations)
    return G_i, eigenvals, converged_eigenvals, H_i


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
                 precomputed_pca=None, federated_qr=False, v=None, gradient=False, epsilon=0.000001):
    """
    Retrieve the first eigenvector
    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of countermensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    """

    iterations = 0  # iteration counter
    converged = False
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
        # log_transmission(filename, "G_i=SC", iterations, i, G_list[i], id+1)
        if V_k is not None:
            Vk_list.append(V_k[start:start + local_data[i].shape[1], :])
        start = start + local_data[i].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)
    G_i_prev = G_i
    while not converged and iterations < maxit:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], 1))  # dummy initialise the H_i matrix
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            # log_transmission(filename, "H_local=CS", iterations, i, H_local ,id+1)
            H_i = H_i + H_local
        # log_transmission(filename, "H_global=SC", iterations, 1, H_i ,id+1)

        for i in range(len(local_data)):
            if gradient:
                G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            else:
                G_list[i] = np.dot(local_data[i].T, H_i)

        gi_norm = 0
        if V_k is None:
            # compute the norm of the eigenvector and done.
            for i in range(len(G_list)):
                local_norm = np.sum(np.square(G_list[i]))
                gi_norm = gi_norm + np.sum(np.square(G_list[i]))
                # log_transmission(filename, "local_norm=CS", iterations, i, local_norm ,id+1)

        elif federated_qr:
            temp = []
            for i in range(len(G_list)):
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
                # log_transmission(filename, "local_dot_prod=CS", iterations, i, np.asarray(sum), id+1)
                local_sums.append(sum)
            local_sums = np.asarray(local_sums)
            local_sums = np.sum(local_sums, axis=0).flatten()  # flatten to remove nesting
            # log_transmission(filename, "global_dot_prod=SC", iterations, 1, local_sums, id+1)

            for i in range(len(G_list)):
                ap = G_list[i]
                for vi in range(Vk_list[i].shape[1]):
                    it = local_sums[vi] * Vk_list[i][:, vi:vi + 1].T
                    it = np.reshape(it, ap.shape)
                    ap = ap - it
                G_list[i] = ap

            for i in range(len(G_list)):
                c_n = np.sum(np.square(G_list[i]))
                # log_transmission(filename, "local_norm=CS", iterations, i, c_n, id+1)
                gi_norm = gi_norm + c_n

        if V_k is None or not federated_qr:
            gi_norm = np.sqrt(gi_norm)
            for i in range(len(G_list)):
                G_list[i] = G_list[i] / gi_norm
                # log_transmission(filename, "ALT_G_i_local=CS", iterations, i, G_list[i], id+1)
                # log_transmission(filename, "ALT_G_i=SC", iterations, i, G_list[i], id+1)

        # log_transmission(filename, "global_norm=SC", iterations, 1, gi_norm, id+1)

        convergedH, deltaH = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        H_i_prev = H_i

        G_i = np.concatenate(G_list, axis=0)
        convergedG, deltaG = sh.eigenvector_convergence_checker(G_i, G_i_prev, tolerance=epsilon)
        G_i_prev = G_i
        log_current_accuracy(u[:, id:id + 1], G_i, [gi_norm], conv=deltaH, current_iteration=iterations,
                             filename=filename, choices=choices, precomputed_pca=precomputed_pca, current_ev=id + 1,
                             gi_delta_obj=deltaG, v=v[:, id:id + 1], H_i=H_i)

    # return a complete eigenvector
    return G_i


def compute_k_eigenvectors(data_list, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None,
                           federated_qr=False, v=None, gradient=True, epsilon=0.000001):
    ug = simulate_guo(data_list, maxit=maxit, filename=filename, u=u, choices=choices, federated_qr=federated_qr, v=v,
                      gradient=gradient, epsilon=epsilon)
    u_all = ug
    for i in range(1, k):
        ug2 = simulate_guo(data_list, maxit=maxit, V_k=u_all, filename=filename, u=u, choices=choices,
                           federated_qr=federated_qr, v=v, gradient=gradient, epsilon=epsilon)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all


def simulate_distributed_horizontal(local, p=8, v=None, choices=None,
                                    precomputed_pca=None, filename=None, maxit=2000, epsilon=0.000001):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of data matrices variable X samples
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[0]
    X = powerit.generate_random_gaussian(d, p, sigma=1)
    X_prev = powerit.generate_random_gaussian(d, p, sigma=1)
    X, R = la.qr(X, mode='economic')
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    while not converged and count < maxit:
        count = count + 1
        # eigenvalues are the column norms of the unnormalised matrix
        X = np.zeros(X.shape[0])
        for l in local:
            #  are in colums.
            d1 = np.dot(l.T, X_prev)
            d1 = np.dot(l, d1)
            X = X + d1
        E = la.norm(X, axis=0)
        X, R = la.qr(X, mode='economic')
        converged, deltas = sh.eigenvector_convergence_checker(X, X_prev, tolerance=epsilon)
        log_current_accuracy(v=v, H_i=X, eigenvals=E, conv=None, current_iteration=count, filename=filename,
                             choices=choices, precomputed_pca=precomputed_pca, u=None, G_i=None)
        X_prev = X
    return X, E, count


####### BENCHMARK RUNNER #######
def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir, epsilon=1e-6, precomputed_pca=None,
                  unequal=False, horizontal=False):
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
            for fedqr, mode in zip([True, False], ['federated_qr', 'central_qr']):

                start = time.monotonic()

                # run power iteration benchmark
                # simultaneous
                grad = False
                grad_name = 'power'
                outdir_gradient = op.join(outdir, 'matrix', str(s), grad_name, mode)
                os.makedirs(outdir_gradient, exist_ok=True)
                filename = create_filename(outdir_gradient, dataset_name + '_' + mode, s, c, k, maxit, start)
                simulate_subspace_iteration(data_list, k, maxit=maxit, u=u, filename=filename, choices=choice,
                                       precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad, epsilon=epsilon)
                end = time.monotonic()
                log_time(logftime, 'qr_scheme' + '_' + mode, end - start, s, c)


                # Run power iteration based benchmark
                # Sequential
                outdir_gradient = op.join(outdir, 'vector', str(s), grad_name, mode)
                os.makedirs(outdir_gradient, exist_ok=True)
                filename = create_filename(outdir_gradient, dataset_name_guo + '_' + mode, s, c, k, maxit, start)

                start = time.monotonic()
                compute_k_eigenvectors(data_list, k=k, maxit=maxit, u=u, filename=filename, choices=choice,
                                       precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad, epsilon=epsilon)
                end = time.monotonic()
                log_time(logftime, 'guo_single' + '_' + mode, end - start, s, c)

                # Run Guo version
                # Sequention
                grad = True
                grad_name = 'gradient'
                outdir_gradient = op.join(outdir, 'vector', str(s), grad_name, mode)
                os.makedirs(outdir_gradient, exist_ok=True)

                filename = create_filename(outdir_gradient, dataset_name_guo + '_' + mode, s, c, k, maxit, start)

                start = time.monotonic()
                compute_k_eigenvectors(data_list, k=k, maxit=maxit, u=u, filename=filename, choices=choice,
                                       precomputed_pca=precomputed_pca, federated_qr=fedqr, v=v, gradient=grad,
                                       epsilon=epsilon)
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
    parser.add_argument('--vert', action='store_true', help='run vertical split test')
    parser.add_argument('--hor', action='store_true', help='run horizontal split test')
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

    print(outdir)
    nr_samples = 0
    nr_features = 0
    if filetype == 'delim-list':
        data_list = []
        for f in path.split(','):
            data, sample_ids, variable_names = si.data_import(f, sep=sep)
            print('tt')
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
                      outdir=vertical, precomputed_pca=precomputed_pca, unequal=unequal)

    # horizontal test
    # transpose data and run vertical test.
    if args.hor:
        data = data.T
        horizontal = op.join(outdir, 'horizontal')
        os.makedirs(horizontal, exist_ok=True)
        the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                      outdir=horizontal, precomputed_pca=precomputed_pca, unequal=unequal, horizontal=True)



    #### TEST CODE ###
    # data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # #
    # #data, sample_ids, variable_names = si.data_import('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t', header=0, rownames=0)
    # #data = si.scale_center_data_columnwise(data)
    # data = coo_matrix.asfptype(data)
    # # args.k = 10
    # # g = gv.standalone(data, k=2)
    # #
    # u, s, v = lsa.svds(data.T,k=10)
    # u = np.flip(u, axis = 1)
    # s = np.flip(s)
    # v = np.flip(v.T, axis=1)
    #
    # data_list, choices = sh.partition_data_vertically(data,3)
    # # ev =compute_k_eigenvectors(data_list, 10, 2000)
    # # ev1 = hybrid_scheme(data_list, 10, 2000, u=u)
    # ev2,ee2,ee3, HI = simulate_subspace_iteration(data_list, 10, 500)
    # #print(co.compute_angles(u, ev))
    # #print(co.compute_angles(u, ev1))
    # angles = co.compute_angles(u, ev2)
    #
    #
    # dsub  = data[0:30000, :]
    # u, s, v = lsa.svds(dsub.T,k=10)
    # u = np.flip(u, axis = 1)
    # s = np.flip(s)
    # v = np.flip(v.T, axis=1)
    #
    # data_list, choices = sh.partition_data_vertically(dsub,3)
    # # ev =compute_k_eigenvectors(data_list, 10, 2000)
    # # ev1 = hybrid_scheme(data_list, 10, 2000, u=u)
    # ev2,ee2,ee3, HI = simulate_subspace_iteration(data_list, 10, 500)
    # #print(co.compute_angles(u, ev))
    # #print(co.compute_angles(u, ev1))
    # angles2 = co.compute_angles(u, ev2)
    #
    #
    # dsub2  = data[0:10000, :]
    # u, s, v = lsa.svds(dsub2.T,k=10)
    # u = np.flip(u, axis = 1)
    # s = np.flip(s)
    # v = np.flip(v.T, axis=1)
    #
    # data_list, choices = sh.partition_data_vertically(dsub2,3)
    # # ev =compute_k_eigenvectors(data_list, 10, 2000)
    # # ev1 = hybrid_scheme(data_list, 10, 2000, u=u)
    # ev2,ee2,ee3, HI = simulate_subspace_iteration(data_list, 10, 500)
    # #print(co.compute_angles(u, ev))
    # #print(co.compute_angles(u, ev1))
    # angles3 = co.compute_angles(u, ev2)
    #
    # # ev = compute_k_eigenvectors(data_list, 10, 2000, federated_qr=True)
    # # ev1 = hybrid_scheme(data_list, 10, 2000, u=u, federated_qr=True)
    # # ev2, ee2, ee3 = simulate_subspace_iteration(data_list, 10, 500, federated_qr=True)
    # # print(co.compute_angles(u, ev))
    # # print(co.compute_angles(u, ev1))
    # # print(co.compute_angles(u, ev2))
    # #
    # #