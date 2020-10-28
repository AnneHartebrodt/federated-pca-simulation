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


####### LOGGING FUNCTIONS #######
def filename(dataset_name, splits, counter, k, maxit, time):
    '''
    Make a file name
    Args:
        dataset_name: Name of the dataset currently running
        splits: The current split
        counter: Current repeat experiment
        k: Number of targeted dimensions
        maxit: maximal iterations

    Returns: A concatenated filename with filestamp

    '''
    fn = dataset_name + '_' + str(splits) + '_' + str(counter) + '_' + str(k) + '_' + str(maxit) + '_' + str(time)
    return fn


def start_logging(outdir, file_ending, dataset_name, maxit, counter, nr_samples, nr_features, k, convergence_eps,
                  splits, time):
    '''

    Args:
        file_ending: File ending of the log file
        dataset_name: Name of the current data set file
        maxit: maximum iterations allowed for this run
        counter: Number of the run, when experiment is repeated
        nr_samples: Number of samples in the dataset (i.e. #patients/indviduals)
        nr_features: NUmber of features in the data set (i.e. #SNPs/genes)
        k: The targeted output dimension

    Returns: The filename of the file without ending

    '''
    fn = filename(dataset_name, splits, counter, k, maxit, time)
    fn = op.join(outdir, fn)
    fn_end = fn + file_ending
    with open(fn_end, 'a') as handle:
        info = '# !data set name\t' + str(dataset_name) + '\n'
        handle.write(info)
        info = '# !data set #features\t' + str(nr_features) + '\n'
        handle.write(info)
        info = '# !data set #samples\t' + str(nr_samples) + '\n'
        handle.write(info)
        info = '# !maximal iterations\t' + str(maxit) + '\n'
        handle.write(info)
        info = '# !target dimensions\t' + str(k) + '\n'
        handle.write(info)
        info = '# !convergence tolerance\t' + str(convergence_eps) + '\n'
        handle.write(info)
        info = '# !run no.\t' + str(counter) + '\n'
        handle.write(info)
        info = '# !data set splits\t' + str(splits) + '\n'
        handle.write(info)
        info = '# !start time\t' + str(time) + '\n'
        handle.write(info)
    return fn


def log_current_accuracy(scipy, G_i, eigenvals, conv, current_iteration, filename, choices, precomputed_pca=None,
                         current_ev=None):
    '''
    Log the current iterations angle to the canonical
    Args:
        scipy: column vector matrix with canonical eigenvectors
        G_i: column vector based matrix with current eigenvector estimation
        current_iteration: iteration index
        filename: output filename prefix > out will be saved to out.angles, and out.cor

    Returns: None

    '''
    if current_ev is not None:
        info_string = str(current_ev) + '\t' + str(current_iteration)
    else:
        info_string = str(current_iteration)

    with open(filename + '.angles', 'a') as handle:
        angles = co.compute_angles(scipy[choices, :], G_i)
        if angles is not None and len(angles) > 0:
            info = cv.collapse_array_to_string(angles, info_string)
            handle.write(info)

    with open(filename + '.cor', 'a') as handle:
        correlations = co.compute_correlations(scipy[choices, :], G_i)
        if correlations is not None and len(correlations) > 0:
            info = cv.collapse_array_to_string(correlations, info_string)
            handle.write(info)

    with open(filename + '.eigenval', 'a') as handle:
        info = cv.collapse_array_to_string(eigenvals, info_string)
        handle.write(info)

    with open(filename + '.conv', 'a') as handle:
        conv = cv.collapse_array_to_string(conv, info_string)
        handle.write(conv)

    if precomputed_pca is not None:
        with open(filename + '.angles_precomp', 'a') as handle:
            angles = co.compute_angles(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(angles, info_string)
            handle.write(info)

        with open(filename + '.cor_precomp', 'a') as handle:
            correlations = co.compute_correlations(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(correlations, info_string)
            handle.write(info)


def log_choices(logfile, filename, choices):
    '''
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    '''
    with open(logfile, 'a') as handle:
        handle.write(cv.collapse_array_to_string(choices, filename))


def log_time(logfile, algorithm, time, split, repeat):
    '''
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    '''
    with open(logfile, 'a') as handle:
        handle.write(algorithm + '\t' + str(split) + '\t' + str(repeat) + '\t' + str(time) + '\n')


def init_benchmark(outdir, dataset_name, maxit, counter,
                   nr_samples, nr_features, k, convergence_eps,
                   splits, timer, precomputed_pca=False, transmission_costs=False):
    filename = start_logging(outdir=outdir, file_ending='.angles', dataset_name=dataset_name, maxit=maxit,
                             counter=counter,
                             nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps,
                             splits=splits, time=timer)
    start_logging(outdir, '.cor', dataset_name=dataset_name, maxit=maxit, counter=counter,
                  nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps,
                  splits=splits, time=timer)
    start_logging(outdir, '.eigenval', dataset_name=dataset_name, maxit=maxit, counter=counter,
                  nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps,
                  splits=splits, time=timer)
    start_logging(outdir, '.conv', dataset_name=dataset_name, maxit=maxit, counter=counter,
                  nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps,
                  splits=splits, time=timer)
    if precomputed_pca is not False:
        start_logging(outdir=outdir, file_ending='.angles_precomp', dataset_name=dataset_name, maxit=maxit,
                      counter=counter, nr_samples=nr_samples, nr_features=nr_features, k=k,
                      convergence_eps=convergence_eps, splits=splits, time=timer)
        start_logging(outdir, '.cor_precomp', dataset_name, maxit, counter, nr_samples, nr_features, k, convergence_eps,
                      splits=splits, time=timer)
    if transmission_costs is not False:
        start_logging(outdir=outdir, file_ending='.transmission', dataset_name=dataset_name, maxit=maxit,
                      counter=counter, nr_samples=nr_samples, nr_features=nr_features, k=k,
                      convergence_eps=convergence_eps, splits=splits, time=timer)
    return filename


####### END LOGGING FUNCTIONS #######

####### MATRIX POWER ITERATION SCHEME #######
def simulate_guo_benchmark(local_data, k, maxit, filename=None, scipy=None, choices=None, precomputed_pca=None, fractev=1.0):
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
        log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
        start = start + local_data[i].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    converged_eigenvals = []
    while not converged and iterations < maxit and len(converged_eigenvals)<k*fractev:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], k))
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            log_transmission(filename, "H_local=CS", iterations, i, H_local)
            H_i = H_i + H_local
        log_transmission(filename, "H_global=SC", iterations, 1, H_i)

        for i in range(len(G_list)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            log_transmission(filename, "Gi_local=CS", iterations, i, G_list[i])

        G_i = np.concatenate(G_list, axis=0)

        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        eigenvals = np.sqrt(eigenvals)

        G_i, R = la.qr(G_i, mode='economic')
        start = 0
        for i in range(len(G_list)):
            G_list[i] = G_i[start:start + local_data[i].shape[1], :]
            log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
            start = start + local_data[i].shape[1]

        converged, conv, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)
        H_i_prev = H_i
        log_current_accuracy(scipy=scipy, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations,
                             filename=filename,
                             choices=choices, precomputed_pca=precomputed_pca)
    G_i = np.concatenate(G_list)
    print(iterations)
    return G_i, eigenvals, converged_eigenvals


####### END MATRIX POWER ITERATION SCHEME #######

def assure_consecutive(arr):
    if len(arr) == 0:
        return -1
    i = 0
    while i < (len(arr) - 1) and arr[i] + 1 == arr[i + 1]:
        i = i + 1
    return i


####### ORIGINAL POWER ITERATION SCHEME #######
def residuals(V, a=None, sums=None):
    '''
    '''
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


def simulate_guo(local_data, maxit, V_k=None, starting_vector=None, filename=None, scipy=None, choices=None,
                 precomputed_pca=None):
    '''
    Retrieve the first eigenvector
    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of countermensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''

    iterations = 0  # iteration counter
    converged = False
    total_len = 0  # total number of samples/incounterviduals in data
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
        log_transmission(filename, "G_i=SC", iterations, i, G_list[i], id+1)
        if V_k is not None:
            Vk_list.append(V_k[start:start + local_data[i].shape[1], :])
        start = start + local_data[i].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not converged and iterations < maxit:

        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], 1))  # dummy initialise the H_i matrix
        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            log_transmission(filename, "H_local=CS", iterations, i, H_local ,id+1)
            H_i = H_i + H_local
        log_transmission(filename, "H_global=SC", iterations, 1, H_i ,id+1)

        for i in range(len(local_data)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]

        gi_norm = 0
        if V_k is None:
            # compute the norm of the eigenvector and done.
            for i in range(len(G_list)):
                local_norm = np.sum(np.square(G_list[i]))
                gi_norm = gi_norm + np.sum(np.square(G_list[i]))
                log_transmission(filename, "local_norm=CS", iterations, i, local_norm ,id+1)

        else:
            local_sums = []
            for i in range(len(G_list)):
                sum = []
                for vi in range(Vk_list[i].shape[1]):
                    dp = np.dot(G_list[i].T, Vk_list[i][:, vi:vi + 1]).flatten()
                    sum.append(dp)
                # cast to numpy array to determine size
                log_transmission(filename, "local_dot_prod=CS", iterations, i, np.asarray(sum), id+1)
                local_sums.append(sum)
            local_sums = np.asarray(local_sums)
            local_sums = np.sum(local_sums, axis=0).flatten()  # flatten to remove nesting
            log_transmission(filename, "global_dot_prod=SC", iterations, 1, local_sums, id+1)

            for i in range(len(G_list)):
                ap = G_list[i]
                for vi in range(Vk_list[i].shape[1]):
                    it = local_sums[vi] * Vk_list[i][:, vi:vi + 1].T
                    it = np.reshape(it, ap.shape)
                    ap = ap - it
                G_list[i] = ap

            for i in range(len(G_list)):
                c_n = np.sum(np.square(G_list[i]))
                log_transmission(filename, "local_norm=CS", iterations, i, c_n, id+1)
                gi_norm = gi_norm + c_n

        gi_norm = np.sqrt(gi_norm)
        for i in range(len(G_list)):
            G_list[i] = G_list[i] / gi_norm
            log_transmission(filename, "ALT_G_i_local=CS", iterations, i, G_list[i], id+1)
            log_transmission(filename, "ALT_G_i=SC", iterations, i, G_list[i], id+1)

        log_transmission(filename, "global_norm=SC", iterations, 1, gi_norm, id+1)

        converged, sum, conv, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)
        H_i_prev = H_i


        G_i = np.concatenate(G_list, axis=0)
        log_current_accuracy(scipy[:,id:id+1], G_i, [gi_norm], conv=delta, current_iteration=iterations,
                             filename=filename, choices=choices, precomputed_pca=precomputed_pca, current_ev=id+1)

    # return a complete eigenvector
    return G_i


def compute_k_eigenvectors(data_list, k, maxit, filename=None, scipy=None, choices=None, precomputed_pca=None):
    ug = simulate_guo(data_list, maxit=maxit, filename=filename, scipy=scipy, choices=choices)
    u_all = ug
    for i in range(1, k):
        ug2 = simulate_guo(data_list, maxit=maxit, V_k=u_all, filename=filename, scipy=scipy, choices=choices)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all


def log_transmission(logfile, log_entry_string, iterations, counter, element, eigenvector=10):
    with open(logfile + '.transmission', 'a') as handle:
        handle.write(log_entry_string + '\t' + str(iterations)
                     + '\t' + str(counter) +'\t' + str(eigenvector)+'\t' + str(sys.getsizeof(element.tobytes()))+'\n')


####### END ORIGINAL POWER ITERATION SCHEME ####


####### HYBRID POWER ITERATION SCHEME #######
def hybrid_scheme(data_list, k, maxit, filename=None, filename2=None, scipy=None, choices=None, precomputed_pca=None):
    #ug, eigenvalues = better_hybrid_scheme(local_data=data_list, k=k, maxit=maxit, filename=filename, choices=choices,
                #                           scipy=scipy, precomputed_pca=precomputed_pca)

    ug, ev, eigenvalues = simulate_guo_benchmark(local_data=data_list, k=k, maxit=maxit, filename=filename, choices=choices, scipy=scipy, precomputed_pca=precomputed_pca, fractev=0.75)
    ug1 = ug.copy()
    restart = assure_consecutive(eigenvalues)
    print(restart)
    for i in range(restart, k - 1):
        next_eigenvec = simulate_guo(data_list, maxit=maxit, V_k=ug1[:, 0:i+1], starting_vector=ug1[:, i+1],
                                     filename=filename2, scipy=scipy, choices=choices, precomputed_pca=None)
        ug = np.concatenate((ug1[:, 0:i+1], next_eigenvec), axis=1)
    return ug


def better_hybrid_scheme(local_data, k, maxit, filename, scipy, choices, precomputed_pca=None):
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
        log_transmission(filename, "G_i=SC", iterations, i, G_list[i])

        start = start + local_data[i].shape[1]

    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)  # dummy init

    G_conv = None  # converged eigenvector matrix
    converged_counter = 0
    annormaly = [0] * k  # init counter for small eigengaps
    global_eigenvalues = []
    global_deltas = []
    while not converged and iterations < maxit and annormaly[converged_counter] < 10:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], k - converged_counter))

        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            log_transmission(filename, "H_local=CS", iterations, i, H_local)
            H_i = H_i + H_local

        log_transmission(filename, "H_global=SC", iterations, 1, H_i)

        converged, sum_of_delta, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)

        for i in range(len(G_list)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]

        G_i = np.concatenate(G_list, axis=0)

        # get the converged eigenvalues
        eigenvals = global_eigenvalues.copy()
        delta = np.append(global_deltas, delta)
        diff = []
        # iterate over non converged eigenvectors
        for col in range(G_i.shape[1]):
            eigenvals.append(np.sqrt(np.linalg.norm(G_i[:, col])))
        # compute eigengaps for all eigenvalues
        for e in range(1, len(eigenvals)):
            diff.append(np.log(np.abs(eigenvals[e - 1] - eigenvals[e])))
        # compute mean and standard deviation over all eigengaps
        mm = np.mean(diff)
        ssd = np.sqrt(np.var(diff))

        # flag eigenvectors that exhibit a very small eigengap
        for d in range(len(diff)):
            if diff[d] < mm - ssd:
                annormaly[d] = annormaly[d] + 1
        #print(annormaly)

        # Orthonormalise based on all eigenvectors
        if G_conv is not None and G_conv.shape[1] > 0:
            G_i = np.concatenate([G_conv, G_i], axis=1)
        G_i, R = la.qr(G_i, mode='economic')

        # current iteration led to some eigenvectors converging
        if len(converged_eigenvals) > 0:
            converged_counter = converged_counter + len(converged_eigenvals)
            # slice H_i and reset annomaly counter
            H_i_prev = H_i[:, len(converged_eigenvals):]
            annormaly = k * [0]
            # update converged eigenvector array
            global_eigenvalues = eigenvals[0:converged_counter]
            global_deltas = delta[0:converged_counter]
        else:
            H_i_prev = H_i

        log_current_accuracy(scipy=scipy, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations,
                             filename=filename, choices=choices, precomputed_pca=precomputed_pca)
        # in any case reslice G_i into converged and
        # unconverged
        G_conv = G_i[:, 0:converged_counter]
        G_i = G_i[:, converged_counter:]

        # redistribute the eigenvector parts
        start = 0
        for i in range(len(local_data)):
            G_list[i] = G_i[start:start + local_data[i].shape[1], :]
            log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
            start = start + local_data[i].shape[1]


    print(iterations)
    if  G_i.shape[1] > 0:
        G_conv = np.concatenate([G_conv, G_i], axis=1)
    return G_conv, global_eigenvalues


####### END HYBRID POWER ITERATION SCHEME #######


####### BENCHMARK RUNNER #######
def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir, convergence_eps=1e-6, precomputed_pca=None):
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
            filename = init_benchmark(outdir=outdir, dataset_name=dataset_name, maxit=maxit, counter=c,
                                      nr_samples=nr_samples, nr_features=nr_features, k=k,
                                      convergence_eps=convergence_eps,
                                      splits=s, timer=timer, transmission_costs=True)

            # split the data
            if not islist:
                data_list, choice = sh.partition_data_vertically(data, s, randomize=True)
            else:
                choice = range(data.shape[1])
            logf = op.join(outdir, 'log_choices.log')
            log_choices(logf, filename, choice)

            start = time.monotonic()
            logftime = op.join(outdir, 'time.log')

            # # simulate the run
            simulate_guo_benchmark(data_list, k, maxit=maxit, scipy=u, filename=filename, choices=choice,
                                   precomputed_pca=precomputed_pca)
            end = time.monotonic()

            log_time(logftime, 'qr_scheme', end - start, s, c)

            filename = init_benchmark(outdir=outdir, dataset_name=dataset_name_guo, maxit=maxit, counter=c,
                                      nr_samples=nr_samples, nr_features=nr_features, k=k,
                                      convergence_eps=convergence_eps,
                                      splits=s, timer=timer, transmission_costs=True)

            start = time.monotonic()
            compute_k_eigenvectors(data_list, k=k, maxit=maxit, scipy=u, filename=filename, choices=choice,
                                            precomputed_pca=precomputed_pca)
            end = time.monotonic()
            log_time(logftime, 'guo_single', end - start, s, c)

            filename = init_benchmark(outdir=outdir, dataset_name=dataset_name + '_hybrid', maxit=maxit, counter=c,
                                      nr_samples=nr_samples, nr_features=nr_features, k=k,
                                      convergence_eps=convergence_eps,
                                      splits=s, timer=timer, transmission_costs=True)
            filename2 = init_benchmark(outdir=outdir, dataset_name=dataset_name + '_hybrid_reit', maxit=maxit,
                                       counter=c,
                                       nr_samples=nr_samples, nr_features=nr_features, k=k,
                                       convergence_eps=convergence_eps,
                                       splits=s, timer=timer, transmission_costs=True)

            start = time.monotonic()
            hybrid_scheme(data_list, k=k, maxit=maxit, scipy=u, filename=filename, filename2=filename2,
                         choices=choice,
                         precomputed_pca=precomputed_pca)
            end = time.monotonic()
            log_time(logftime, 'hybrid_scheme', end - start, s, c)


####### BENCHMARK RUNNER #######

if __name__ == '__main__':
   #  data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
   #  data = coo_matrix.asfptype(data)
   #  # args.k = 10
   #  # g = gv.standalone(data, k=2)
   #  #
   #  u, s, v = lsa.svds(data.T,k=10)
   #  u = np.flip(u, axis = 1)
   #  s = np.flip(s)
   #  v = np.flip(v.T, axis=1)
   #
   #  data_list, choices = sh.partition_data_vertically(data,2)
   #  ev = hybrud(data_list, 10, 2000, scipy=u)
   #  print(co.compute_angles(u, ev))
   # # print(co.angle(u[:,7], ev[:,9]))




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
    parser.add_argument('--lvg', action='store_true', help='spreadsheet separator, default tab')
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

    s = args.s
    splits = s.strip().split(',')
    splits = np.int8(splits)

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

    the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits,
                  outdir=outdir, precomputed_pca=precomputed_pca)

    # produce k-1 eigenvectors

    if filetype == 'delim' and args.orthovector is not None:
        # data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw','train')
        data_list, choices = sh.partition_data_vertically(data, 2)
        ug, ev = runner.simulate_guo(data_list, k, maxit=maxit)

        ortho = []
        aps = []
        for i in range(100):
            ap = gv.get_initial_eigenvector_k(ug)
            aps.append(ap)
            loca = []
            for v in range(ug.shape[1]):
                loca.append(co.angle(ug[:, v], ap.T))
            ortho.append(loca)

        ap_angles = []
        for a in range(len(aps) - 1):

            for a1 in range(a + 1, len(aps)):
                print(co.angle(ap[a], ap[a1].T))
                ap_angles.append(co.angle(ap[a], ap[a1]))

        ortho = np.asarray(ortho)

        pd.DataFrame(ortho).to_csv(args.orthovector, header=False)

