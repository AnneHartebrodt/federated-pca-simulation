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
import os.path as op
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


def simulate_guo_benchmark(local_data, k, maxit, filename, scipy, choices, precomputed_pca=None):
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
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    # send parts to local sites
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

        ra, conv = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
        log_current_accuracy(scipy, G_i, eigenvals, conv, current_iteration=iterations, filename=filename, choices=choices, precomputed_pca=precomputed_pca)

    G_i = np.concatenate(G_list)
    return G_i, eigenvals

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

def start_logging(outdir, file_ending, dataset_name, maxit, counter, nr_samples, nr_features, k, convergence_eps, splits, time):
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
    fn_end = fn+file_ending
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
        info = '# !start time\t' + str(time)+'\n'
        handle.write(info)
    return fn


def log_current_accuracy(scipy, G_i, eigenvals, conv, current_iteration, filename, choices, precomputed_pca=None):
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

    with open(filename+'.conv', 'a') as handle:
        info = str(current_iteration)+'\t'+str(conv)+'\n'
        handle.write(info)

    if precomputed_pca is not None:
        with open(filename + '.angles_precomp', 'a') as handle:
            angles = co.compute_angles(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(angles, str(current_iteration))
            handle.write(info)

        with open(filename + '.cor_precomp', 'a') as handle:
            correlations = co.compute_correlations(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(correlations, str(current_iteration))
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

def the_epic_loop(data, dataset_name, maxit, nr_repeats, k, splits, outdir, precomputed_pca=None):
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
    #g = gv.standalone(data, k)
    if isinstance(data, list):
        data_list = data
        data = np.concatenate(data, axis=1)
        splits =[1] # data is already split, only counter experiments need to be run.

    u, s, v = lsa.svds(data.T, k=k)
    u = np.flip(u, axis=1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    convergence_eps = 0.00000001

    for c in range(nr_repeats):
        for s in splits:
            # filename will be the same for angle log file and correlation log file
            timer = time.monotonic()
            filename = start_logging(outdir=outdir, file_ending='.angles', dataset_name=dataset_name, maxit=maxit, counter=c, nr_samples=nr_samples, nr_features=nr_features, k=k, convergence_eps=convergence_eps, splits=s, time = timer)
            start_logging(outdir, '.cor', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps, splits=s, time = timer)
            start_logging(outdir, '.eigenval', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps,
                          splits=s,  time = timer)
            start_logging(outdir, '.conv', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps,
                          splits=s, time = timer)
            if precomputed_pca is not None:
                filename = start_logging(outdir=outdir, file_ending='.angles_precomp', dataset_name=dataset_name, maxit=maxit,
                                         counter=c, nr_samples=nr_samples, nr_features=nr_features, k=k,
                                         convergence_eps=convergence_eps, splits=s)
                start_logging(outdir, '.cor_precomp', dataset_name, maxit, c, nr_samples, nr_features, k, convergence_eps,
                              splits=s)


            # split the data
            if not isinstance(data, list):
                data_list, choice = sh.partition_data_vertically(data, s, randomize=True)

            logf = op.join(outdir, 'log_choices.log')
            log_choices(logf, filename, choice)

            # simulate the run
            simulate_guo_benchmark(data_list, k + 2, maxit=maxit, scipy=u, filename=filename, choices=choice, precomputed_pca=precomputed_pca)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; default tab separated')
    parser.add_argument('--filetype', metavar='filetype', type=str, help='Type of the dataset')
    parser.add_argument('--sep', metavar='sep', type=str, help='spreadsheet separator, default tab', default='\t')
    parser.add_argument('--variance', action='store_true', help='spreadsheet separator, default tab')
    parser.add_argument('--center', action='store_true',help='center data')
    parser.add_argument('-o', metavar='outfile', type=str, help='output directory')
    parser.add_argument('-r', metavar='repeats', type=int, default=20, help = 'Number of times to repeat experiment')
    parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
    parser.add_argument('-s', metavar='sites', default='2,3,5,10', type=str, help='comma separated list of number of sites to simulate, parsed as string')
    parser.add_argument('-i', metavar='iteration', default=2000, type=int, help='Maximum number of iterations')
    parser.add_argument('--header', metavar='iteration', default=None, type=int, help='header lines')
    parser.add_argument('--rownames', metavar='iteration', default=None, type=int, help='rownames')
    parser.add_argument('--names', metavar='iteration', default=None, type=str, help='names')
    parser.add_argument('--compare_pca', metavar='compare', default=None, type=str, help='filename of precomputed pca to be compared to')
    parser.add_argument('--orthovector', metavar='compare', default=None, type=str,
                        help='filename of orthogonal file')
    args = parser.parse_args()

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
    #orthovector = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/angles_ortho.tsv'

    print(outdir)
    nr_samples = 0
    nr_features = 0
    if filetype =='delim-list':
        data_list = []
        for f in path.split(','):
            data, sample_ids, variable_names = si.data_import(f, sep=sep )
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
        bim = path+'.bim'
        traw = path+'.traw'
        traw_nosex = gi.remove_non_autosomes(bim, traw)
        data = gi.read_scale_write(infile=traw_nosex, outfile=path+'.traw.scaled', maf=0.01)
        nr_samples = data.shape[0]
        nr_features = data.shape[1]
    else:
        raise Exception("Filetype not supported")

    if args.compare_pca is not None:
        precomputed_pca = pd.read_table(args.compare_pca, header=0, sep='\t')
        precomputed_pca = precomputed_pca.values
    else:
        precomputed_pca=None


    the_epic_loop(data=data, dataset_name=dataset_name, maxit=maxit, nr_repeats=nr_repeats, k=k, splits=splits, outdir=outdir, precomputed_pca=precomputed_pca)


    # produce k-1 eigenvectors

    if filetype == 'delim' and args.orthovector is not None:
        #data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw','train')
        data_list, choices = sh.partition_data_vertically(data, 2)
        ug, ev = runner.simulate_guo(data_list, 12, maxit=200)


        ortho = []
        aps = []
        for i in range(100):
            ap = gv.get_initial_eigenvector_k(ug)
            aps.append(ap)
            loca = []
            for v in range(ug.shape[1]):
                loca.append(co.angle(ug[:,v], ap.T))
            ortho.append(loca)

        ap_angles = []
        for a in range(len(aps)-1):
            print(a)
            for a1 in range(a+1, len(aps)):
                print(co.angle(ap[a], ap[a1].T))
                ap_angles.append(co.angle(ap[a], ap[a1]))

        ortho = np.asarray(ortho)

        pd.DataFrame(ortho).to_csv(args.orthovector, header=False)


