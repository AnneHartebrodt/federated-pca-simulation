'''
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

'''

import python.PCA.guo_vertical as gv
import python.PCA.shared_functions as sh
import scipy.linalg as la
#import import_export.easy_import as easy
import scipy.sparse.linalg as lsa
import argparse as ap
import pandas as pd
import os.path as path
import python.PCA.comparison as co
import python.PCA.convenience as cv
import time as time
import python.import_export.mnist_import as imnist
import python.PCA.guo_vertical_runner as runner

from scipy.sparse import coo_matrix
import os
import gzip
import numpy as np
import python.import_export.spreadsheet_import as si
import python.import_export.gwas_import as gi
import python.import_export.mnist_import as mi

# Copied from fashion mnist


def simulate_guo_benchmark(local_data, k, maxit, filename, scipy, choices):
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
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    while not ra and iterations<maxit:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], k))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            #G_i, Q = la.qr(G_i, mode='economic')
            G_list_n.append(G_i)
        start = 0
        G_i = np.concatenate(G_list_n, axis=0)
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        eigenvals = np.sqrt(eigenvals)
        G_i, R = la.qr(G_i, mode='economic')
        G_list = []
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]

        ra = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
        log_current_accuracy(scipy, G_i, eigenvals, current_iteration=iterations, filename=filename, choices=choices)

    G_i = np.concatenate(G_list)
    return G_i, eigenvals

def filename(dataset_name, splits, counter, k, maxit):
    fn = dataset_name + '_' + str(splits) + '_' + str(counter) + '_' + str(k) + '_' + str(maxit)
    return fn

def start_logging(outdir, file_ending, dataset_name, maxit, counter, nr_samples, nr_features, k, convergence_eps, splits):
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
    fn = filename(dataset_name, splits, counter, k, maxit)
    fn = path.join(outdir, fn)
    fn_end = fn+file_ending
    fn_end = path.join(outdir, fn_end)
    with open(fn_end, 'a') as handle:
        info = '# !data set name\t' + str(dataset_name)+'\n'
        handle.write(info)
        info = '# !data set #features\t' + str(nr_features)+'\n'
        handle.write(info)
        info = '# !data set #samples\t' + str(nr_samples)+'\n'
        handle.write(info)
        info = '# !maximal iterations\t'+str(maxit)+'\n'
        handle.write(info)
        info = '# !target dimensions\t' + str(k)+'\n'
        handle.write(info)
        info = '# !convergence tolerance\t' + str(convergence_eps)+'\n'
        handle.write(info)
        info = '# !run no.\t'+str(counter)+'\n'
        handle.write(info)
        info = '# !data set splits\t' + str(splits) + '\n'
        handle.write(info)
        info = '# !start time\t' + str(time.monotonic())+'\n'
        handle.write(info)
    return fn


def log_current_accuracy(scipy, G_i, eigenvals, current_iteration, filename, choices):
    '''
    Log the current iterations angle to the canonical
    Args:
        scipy: column vector matrix with canonical eigenvectors
        G_i: column vector based matrix with current eigenvector estimation
        current_iteration: iteration index
        filename: output filename prefix > out will be saved to out.angles, and out.cor

    Returns: None

    '''
    with open(filename+'.angles', 'a') as handle:
        angles = co.compute_angles(scipy[choices, :], G_i)
        info = cv.collapse_array_to_string(angles, str(current_iteration))
        handle.write(info)

    with open(filename+'.cor', 'a') as handle:
        correlations = co.compute_correlations(scipy[choices, :], G_i)
        info = cv.collapse_array_to_string(correlations, str(current_iteration))
        handle.write(info)

    with open(filename+'.eigenval', 'a') as handle:
        info = cv.collapse_array_to_string(eigenvals, str(current_iteration))
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

def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir):
    '''
    run the simulation of a federated run of vertical power iteration
    Args:
        data:
        dataset_name:
        maxit:
        nr_repeats:
        k:
        splits:
        outdir:

    Returns:

    '''
    #g = gv.standalone(data, k)

    u, s, v = lsa.svds(data.T, k=k)
    u = np.flip(u, axis=1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    convergence_eps = 0.00000001
    print(splits)
    for c in range(nr_repeats):
        for s in splits:
            print(s)
            # filename will be the same for angle log file and correlation log file
            filename = start_logging(outdir=outdir, file_ending='.angles', dataset_name=dataset_name, maxit=maxit, counter=c, nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps, splits=s)
            start_logging(outdir, '.cor', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps, splits=s)
            start_logging(outdir, '.eigenval', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps,
                          splits=s)

            # split the data
            data_list, choice = sh.partition_data_vertically(data, s, randomize=True)
            logf = path.join(outdir, 'log_choices.log')
            log_choices(logf, filename, choice)

            # simulate the run
            simulate_guo_benchmark(data_list, k + 2, maxit=maxit, scipy=u, filename=filename, choices=choice)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; default tab separated')
    parser.add_argument('--filetype', metavar='filetype', type=str, help='Type of the dataset')
    parser.add_argument('--sep', metavar='sep', type=str, help='spreadsheet separator, default tab', default='\t')
    parser.add_argument('-o', metavar='outfile', type=str, help='output directory')
    parser.add_argument('-r', metavar='repeats', type=int, default=20, help = 'Number of times to repeat experiment')
    parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
    parser.add_argument('-s', metavar='sites', default='2,3,5,10', type=str, help='comma separated list of number of sites to simulate, parsed as string')
    parser.add_argument('-i', metavar='iteration', default=2000, type=int, help='Maximum number of iterations')
    args = parser.parse_args()

    # import scaled SNP file
    #data = easy.easy_import(args.f, header=None, rownames=None, center=False, scale_var=False,sep='\t')
    #data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    #data = coo_matrix.asfptype(data)

    path = args.file
    filetype = args.filetype
    sep = args.sep
    k = args.k
    dataset_name = os.path.basename(args.file)
    s = args.s
    splits = s.strip().split(',')
    maxit = args.i
    nr_repeats = args.r
    outdir = args.o

    if filetype == 'delim':
        data = si.data_import(path, sep=sep)
    elif filetype == 'mnist':
        data, test_lables = mi.load_mnist(path, 'train')
    elif filetype == 'gwas':
        data = gi.import_bed(path, True)
    else:
        raise Exception("Filetype not supported")

    nr_samples = data.shape[0]
    nr_features = data.shape[1]

    the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits, outdir=outdir)





