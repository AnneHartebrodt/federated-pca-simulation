import os.path as path
import time as time
import argparse as ap
import numpy as np
import pandas as pd

import comparison as co
import convenience as cv
import import_export.easy_import as easy
import outlier_removal as orm
import power_iteration_runner as power_runner
import proxy_covariance as dpca
import proxy_covariance_runner as proxy_runner


def unqeal_split_power_iteration(data, interval_end, p=10, tolerance=0.000001):
    # Shuffling data to create random samples
    np.random.shuffle(data)
    start = 0
    localdata = []
    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[start:end, :]
        # calculate covariance matrix
        noisy_cov = dpca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, noise=False)
        start = int(interval_end[i])
        localdata.append(noisy_cov)
    W, E, count = power_runner.simulate_distributed(localdata, p)
    return W, E, count


def local_outlier_removal(data_sub, pca, study_id=None, outfile=None):
    print('Local outlier removal')
    outliers = orm.outlier_removal_mad(pca, 6, 3)
    if outfile is not None and study_id is not None:
        with open(path.join(outfile, 'removed_outliers.tsv'), 'a+') as handle:
            handle.write(cv.collapse_array_to_string(outliers, study_id))
    if len(outliers) != 0:
        data_sub = np.delete(data_sub, outliers, 0)
    return data_sub


def unqeal_split_proxy_covariance(data, interval_end, ndims=100, exp_var=0.5, mult_dims_ret=1,
                                  weights=None, balacan=False, unweighted=False, remove_local_outliers=False):
    """
    This function simulates a multisite PCA with each site having
    varying number of samples.
    :param data:
    :param sites: number of sites to split the data into
    :return:
    """

    print('Shuffling data')
    np.random.shuffle(data)

    Ac = []
    variance_explained = []
    start = 0
    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[start:end, :]
        if remove_local_outliers:
            data_sub = local_outlier_removal(data_sub, ndims)
        noisy_cov = dpca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, noise=False)
        start = int(interval_end[i])
        # return the local PCA with the maximal number of dimensions required.
        A, vex = dpca.perform_SVD(noisy_cov, ndims=ndims, mult_dims_returned=max(mult_dims_ret), var_exp=exp_var)
        Ac.append(A)
        variance_explained.append(vex)

    eigenvectors = {'unweighted': [], 'weighted': [], 'balcan': []}
    eigenvalues = {'unweighted': [], 'weighted': [], 'balcan': []}
    vex_max = max(variance_explained)
    print('Aggregate local PCAs')
    for mm in mult_dims_ret:
        mm = int(np.ceil(vex_max * mm))
        if unweighted:
            W, X = dpca.aggregate_partial_SVDs(Ac, ndim=ndims, intermediate_dims=mm)
            eigenvectors['unweighted'].append(W)
            eigenvalues['unweighted'].append(X)
        if weights is not None:
            W, X = dpca.aggregate_partial_SVDs(Ac, ndim=ndims, intermediate_dims=mm, weights=weights)
            eigenvectors['weighted'].append(W)
            eigenvalues['weighted'].append(X)
        if balacan:
            W, X = dpca.aggregate_partial_SVDs_balacan(Ac, ndim=ndims, intermediate_dims=mm)
            eigenvectors['balcan'].append(W)
            eigenvalues['balcan'].append(X)
    return eigenvectors, eigenvalues


def cluster_split(data, clusterfile, ndims=100, p=20, header_clu=None, sep='\t', exp_var=0.5,
                  remove_local_outliers=False,
                  proxy_cov=False, power_it=False):
    """
    This function simulates a multisite PCA with each site having
    varying number of samples.
    :param data:
    :param sites: number of sites to split the data into
    :return:
    """

    Ac = []
    weights = []
    clusters = pd.read_csv(filepath_or_buffer=clusterfile, header=header_clu, sep=sep)
    covs = []
    for i in range(max(clusters.iloc[:, 1])):
        index = clusters.iloc[:, 1] == (i + 1)
        index = clusters[index].index
        # slice matrix
        data_sub = data[index, :]
        if remove_local_outliers:
            data_sub = local_outlier_removal(data_sub, ndims)

        noisy_cov = dpca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, noise=False)
        if power_it:
            covs.append(noisy_cov)
        if proxy_cov:
            A, vex = dpca.perform_SVD(noisy_cov, ndims=ndims, mult_dims_returned=ndims, var_exp=exp_var)
            Ac.append(A)
            weights.append(len(index) * 1.0 / clusters.shape[0])

    eigenvectors = {'powerit': [], 'weighted': []}
    eigenvalues = {'powerit': [], 'weighted': []}

    eigenvectors['weighted'], eigenvalues['weighted'] = dpca.aggregate_partial_SVDs(Ac, ndim=ndims, weights=weights)
    eigenvectors['powerit'], eigenvalues['powerit'], nr_iter = power_runner.simulate_distributed(covs, p)

    return eigenvectors, eigenvalues, nr_iter


def time_logger(task, start, filename=None):
    current = time.monotonic()
    elapsed = current - start
    print(task + '\t' + ' ...time elapsed: ' + str(elapsed))
    if filename is not None:
        with open(path.join(filename, 'time_log.tsv'), 'a+') as handle:
            handle.write(task + '\t' + str(elapsed) + '\n')
    return time.monotonic()


def single_site(data, study_id, dims=100, p=20, outfile=''):
    # single site pca before outlier removal
    start = time.monotonic()
    pca, W0, E0 = proxy_runner.run_standalone(data, dims=dims)
    start = time_logger('Single site PCA', start, outfile)
    outlier_free = local_outlier_removal(data_sub=data, pca=pca)
    start = time_logger('Outlier removal', start, outfile)
    # single after outlier removal
    pca, W1, E1 = proxy_runner.run_standalone(outlier_free, dims=dims)
    start = time_logger('Single site PCA', start, outfile)
    # subspace iteration before outlier removal
    pca, W2, E2, nriter, params = power_runner.run_standalone(data, p=p)
    start = time_logger('Subspace iteration', start, outfile)
    # single site power iteration with deflation
    W3, E3, nr_iter = power_runner.simulate_deflation(data, dims=p)
    start = time_logger('Power iteration', start, outfile)
    # save number of iterations
    with open(path.join(outfile, 'nr_iter_powit.tsv'), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(nr_iter, study_id))
    # number of variables to achieve certain explanatory power
    with open(path.join(outfile, 'nr_vars_explain_aor.tsv'), 'a+') as handle:
        for var in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            handle.write(str(var) + '\t' + str(dpca.variance_explained(E1, perc=var)) + '\n')
    # save number of required iterations until convergences
    with open(path.join(outfile, 'nr_iter_qr.tsv'), 'a+') as handle:
        handle.write(study_id + '\t' + str(nriter) + '\n')

    start = time_logger('Writing files', start, outfile)
    dw = {'single_site_bor': W0, 'single_site_aor': W1, 'single_site_subspace': W2, 'single_site_deflation': W3}
    de = {'single_site_bor': E0, 'single_site_aor': E1, 'single_site_subspace': E2, 'single_site_deflation': E3}
    return dw, de


def write_eigenvectors_single_site(de, reported_angles, outfile):
    for key in de.keys():
        with open(path.join(outfile, key + '_eigenvalues.tsv'), 'a+') as handle:
            handle.write(cv.collapse_array_to_string(de[key][0:reported_angles], study_id))


def write_angles_single_site(dw, reported_angles, outfile):
    keys = list(dw.keys())
    for k1 in range(len(dw.keys())):
        for k2 in range(k1 + 1, len(dw.keys())):
            co.compute_save_angles(dw[keys[k1]], dw[keys[k2]], study_id=study_id,
                                   filename=keys[k1] + '_' + keys[k2] + '_angles.tsv',
                                   outfile=outfile, reported_angles=reported_angles)


def write_single_site(dw, de, reported_angles, outfile):
    write_angles_single_site(dw, reported_angles, outfile=outfile)
    write_eigenvectors_single_site(de, reported_angles, outfile=outfile)


def run_and_compare_unequal(data, outfile, dims=100, p=-1, clusterfile=None, cluster_sep='\t', study_id='',reported_angles=20, exp_var=0.5, mult_dims_ret=[0.5, 1, 2], balcan=False, unweighted = False, weighted = False):
    n = data.shape[0]
    dims = min(dims, n)
    interval_end, perc = make_test_intervals(n)

    dw, de = single_site(data, study_id, p=p, dims=dims, outfile=outfile)
    write_single_site(dw, de, reported_angles, outfile=outfile)

    for ar in range(len(interval_end)):
        for i in range(10):
            print('Current split ' + str(interval_end[ar]))
            start = time.monotonic()
            if unweighted or balcan or weighted:
                if not weighted:
                    perc = None
                eigenvectors_prox, eigenvalues_prox = unqeal_split_proxy_covariance(data, interval_end[ar], ndims=dims,mult_dims_ret=mult_dims_ret,exp_var=exp_var, weights=perc[ar],balacan=balcan, unweighted=unweighted)

                start = time_logger('Unequal split proxy', start, outfile)
                write_results_prox(eigenvectors_prox=eigenvectors_prox, eigenvalues_prox=eigenvalues_prox,
                               reference=dw['single_site_bor'], mult_dims_ret=mult_dims_ret,
                               reported_angles=reported_angles, study_id=study_id, it=i, outfile=outfile)
            if p!=-1:
                eigenvectors_pit, eigenvalues_pit, count_pit = unqeal_split_power_iteration(data, interval_end[ar], p)
                start = time_logger('Unequal split subspace iteration', start, outfile)
                write_results(eigenvectors_pit=eigenvectors_pit, reference=dw['single_site_subspace'],
                          eigenvalues_pit=eigenvalues_pit, study_id=study_id, reported_angles=reported_angles,
                          it=i, file_id='single_site_subspace_dpit_', outfile=outfile)

            # create and write metadata
            meta = [i] + [len(interval_end[ar])] + interval_end[ar]
            with open(path.join(outfile, 'meta_splits.tsv'), 'a+') as handle:
                handle.write(cv.collapse_array_to_string(meta, study_id))

            if clusterfile is not None:
                try:
                    eigenvectors_ueq, eigenvalue_ueq, nr_iter = cluster_split(data, clusterfile=clusterfile,
                                                                          sep=cluster_sep,
                                                                          proxy_cov=True, power_it=True)
                    start = time_logger('Unequal split preclustered', start, outfile)
                    write_results(eigenvectors_pit=eigenvectors_ueq['weighted'], reference=dw['single_site_bor'],
                              eigenvalues_pit=eigenvalue_ueq['weighted'],
                              study_id=study_id,
                              reported_angles=reported_angles, it=i, file_id='single_site_bor_cluster_weighted_',
                              outfile=outfile)
                    write_results(eigenvectors_pit=eigenvectors_ueq['powerit'], reference=dw['single_site_subspace'],
                              eigenvalues_pit=eigenvalue_ueq['powerit'],
                              study_id=study_id, reported_angles=reported_angles, it=i, nr_it=nr_iter,
                              file_id='single_site_bor_cluster_weighted_', outfile=outfile)
                except FileNotFoundError:
                    print('File does not exist. Are your sure there is a preclustered file for this dataset?')




def write_results_prox(outfile, eigenvectors_prox, eigenvalues_prox, reference, mult_dims_ret, reported_angles,
                       study_id, it):
    for key in eigenvectors_prox.keys():
        # list of eigenvector matrices of length mult_dims_ret
        for w in range(len(eigenvectors_prox[key])):
            angles = co.compute_angles(eigenvectors_prox[key][w], reference, reported_angles=reported_angles)
            with open(path.join(outfile, 'prox_angles_unequal_splits_' + key + '_' + str(mult_dims_ret[w]) + '.tsv'),
                      'a+') as handle:
                handle.write(cv.collapse_array_to_string(angles, study_id=study_id))
            with open(path.join(outfile, 'prox_eigenvalues_' + key + '_' + str(mult_dims_ret[w]) + '.tsv'),
                      'a+') as handle:
                handle.write(cv.collapse_array_to_string(eigenvalues_prox[key][0:reported_angles], str(it)))


def write_results(outfile, eigenvectors_pit, reference, eigenvalues_pit, study_id, reported_angles, it, nr_it=-1,
                  file_id=''):
    co.compute_save_angles(eigenvectors_pit, reference, study_id=study_id,
                           filename=file_id + 'angles_unequal_splits.tsv',
                           outfile=outfile, reported_angles=reported_angles)
    with open(path.join(outfile, file_id + 'eigenvalues.tsv'), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(eigenvalues_pit[0:reported_angles], str(it)))
    if nr_it != -1:
        with open(path.join(outfile, file_id + 'iterations_until_convergence.tsv'), 'a+') as handle:
            handle.write(study_id + '\t' + str(it) + '\t' + str(nr_it))


def make_test_intervals(n):
    # hardcoded for now
    # more extreme cases maybe later
    unequal_splits = [[0.1, 0.9], [0.3, 0.7], [0.5, 0.5], [0.2, 0.2, 0.2, 0.2, 0.2], [0.1, 0.1, 0.2, 0.2, 0.4],
                      [0.1, 0.1, 0.1, 0.1, 0.6], [0.2375, 0.2375, 0.2375, 0.2375, 0.05]]
    interval_end = list()
    sum = 0
    for i in unequal_splits:
        inter = list()
        for k in i:
            sum = sum + k
            inter.append(np.ceil(n * sum))
        sum = 0
        interval_end.append(inter)
    return interval_end, unequal_splits


def parse_array(value_str):
    values_str = value_str.split(',')
    values_int = []
    for v in values_str:
        values_int.append(float(v))
    return values_int


if __name__ == "__main__":
    print('run split script')

    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; default tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output directory')
    parser.add_argument('-c', metavar='clusterfile', type=str, help='clustersplit file: tab separated by default')
    parser.add_argument('-n', metavar='header', type=int, help='has data column names?', default=None)
    parser.add_argument('-r', metavar='rownames', type=int, help='has data row names?', default = None)
    parser.add_argument('-s', metavar='sep', type=str, help='field delimiter input file', default = '\t')
    parser.add_argument('-S', metavar='cluster_sep', type=str, help='field delimiter cluster file', default='\t')
    parser.add_argument('-v', metavar='explained_var', type=float, help='explained variance')
    parser.add_argument('-m', metavar='mult_dims_ret', type=str, help='comma separated list of intermediate dimensions',default='1')
    parser.add_argument('-d', metavar='dims', type=int, help='intermediate dimensions single site/proxy pca', default=100)
    parser.add_argument('-p', metavar='powdim', type=int, help='# eigenvectors poweriteration / reported dimensions', default=20)
    parser.add_argument('--center', action='store_true')
    parser.add_argument('--log2', action='store_true')
    parser.add_argument('--balcan', action='store_true')
    parser.add_argument('--weighted', action='store_true')
    parser.add_argument('--unweighted', action='store_true')

    args = parser.parse_args()

    inputfile = args.f
    outfile = args.o
    clusterfile = args.c
    header = args.n
    rownames = args.r
    sep = args.s
    cluster_sep = args.S
    exp_var = args.v
    mult_dims_ret = args.m
    dims  = args.d
    p = args.p
    center = args.center
    log = args.log2
    weighted = args.weighted
    balcan = args.balcan
    unweighted = args.unweighted



    # cluster_sep = '\t'
    # exp_var = 0.5
    # sep = '\t'
    # mult_dims_ret = '0.75,5,10'
    # dims = 100
    # header = 0
    # center = True
    # log = True
    # rownames = None
    # p = 10
    # inputfile = '/home/anne/Documents/featurecloud/data/tcga/data_clean/CPTAC-2/coding_trunc.tsv'
    # outfile = '/home/anne/Documents/featurecloud/results/sandbox2/'
    # clusterfile = '/home/anne/Documents/featurecloud/results/pca_plots/cluster/CPTAC-2_10_clusters.tsv'
    # weighted = True
    # balcan = True
    # unweighted = True

    mult_dims_ret = parse_array(mult_dims_ret)
    summaryfile = cv.make_eigenvector_path(outfile, path.join(path.basename(path.dirname(inputfile)), str(exp_var)))
    study_id = path.basename(path.dirname(inputfile))

    st = time.monotonic()
    # import data, center (can be done globally in the distributed case), log2 transform
    data = easy.easy_import(inputfile, header=header, rownames=rownames, center=center, log=log, sep=sep)

    run_and_compare_unequal(data, summaryfile, reported_angles=p, study_id=study_id, exp_var=exp_var,mult_dims_ret=mult_dims_ret, clusterfile=clusterfile, cluster_sep=cluster_sep, dims=dims,p=p, weighted=weighted, unweighted=unweighted, balcan=balcan)
    time_logger("Total time", st, filename=outfile)