import os.path as path
import comparison as co
import numpy as np

import proxy_covariance_runner as runner
import power_iteration_runner as powerit
import convenience as cv
import outlier_removal as orm
import power_iteration as power
from import_export.import_data import CustomDataImporter
import proxy_covariance as dpca
import import_export.easy_import as easy


def unqeal_split(data, interval_end, p=10, tolerance=0.000001):
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
    W, E, count = powerit.simulate_distributed(localdata, p)
    return W, E, count


def run_and_compare_unequal(datafile, outfile=None, dims=10, header=0, rownames=0, center=True,
                            scale_var=False, scale01=False, scale_unit=False, transpose=False, sep='\t',
                            reported_angles=20, exp_var=0.5, log=True):
    # ID is the folder name
    study_id = path.basename(path.dirname(datafile))
    data = easy.easy_import(datafile, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,
                            scale01=scale01, scale_unit=scale_unit, log=log, outfile=outfile,
                            transpose=transpose)

    n = data.shape[0]
    dims = min(dims, n)

    print('Standalone PCA before outlier removal')
    # Standalone PCA before outlier removal
    pca, W0, E0 = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,
                                        center=False, scale_var=False, scale01=False, scale_unit=False,
                                        sep=sep, log=False, exp_var=exp_var)

    # Logging outliers
    outliers = orm.outlier_removal_mad(pca, 6, 3)
    with open(outfile + '/removed_outliers.tsv', 'a+') as handle:
        handle.write(cv.collapse_array_to_string(outliers, study_id))

    print('Standalone PCA after outlier removal')
    # Standalone PCA after outlier removal
    pca, W0, E0 = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,
                                        center=False, scale_var=False, scale01=False, scale_unit=False,
                                        sep=sep,
                                        drop_samples=outliers, log=False, exp_var=exp_var)

    # Standalone subspace iterateion
    pca, W1, E1, nriter, params = powerit.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,
                                                 center=False, scale_var=False, scale01=False,
                                                 scale_unit=False, sep=sep, log=False, exp_var=exp_var)
    with open(path.join(outfile, 'nr_iter_qr.tsv'), 'a+') as handle:
        handle.write(study_id+'\t'+str(nriter)+'\n')
    # standalone power iteration with deflation
    W2, E2, nr_iter = powerit.simulate_deflation(data, outfile, dims=dims, header=header, rownames=rownames,
                                        center=False, scale_var=False, scale01=False, scale_unit=False,
                                        sep=sep, log=False, exp_var=exp_var)
    # save number of iterations
    with open(path.join(outfile, 'nr_iter_powit.tsv'), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(nr_iter, study_id))

    co.compute_save_angles(W1, W0, study_id=study_id, outfile=outfile, filename='angles_unequal_splits_qr_reg.tsv',
                           reported_angles=reported_angles)
    co.compute_save_angles(W1, W2, study_id=study_id, outfile=outfile,
                           filename='angles_unequal_splits_qr_deflation.tsv',
                           reported_angles=reported_angles)
    co.compute_save_angles(W2, W0, study_id=study_id, outfile=outfile,
                           filename='angles_unequal_splits_reg_deflation.tsv',
                           reported_angles=reported_angles)

    interval_end, perc = make_test_intervals(n)

    for ar in range(len(interval_end)):
        for i in range(1):
            print('Current split' + str(interval_end[ar]))
            W3, E3, count = unqeal_split(data, interval_end[ar], p=dims)

            co.compute_save_angles(W3, W0, study_id=study_id, filename='angles_unequal_splits_split_reg.tsv',
                                   outfile=outfile,
                                   reported_angles=reported_angles)
            co.compute_save_angles(W1, W3, study_id=study_id, outfile=outfile,
                                   filename='angles_unequal_splits_qr_split.tsv',
                                   reported_angles=reported_angles)

            # write results

            meta = [count]+[len(interval_end[ar])] + interval_end[ar]
            with open(path.join(outfile, 'meta_splits.tsv'), 'a+') as handle:
                # #run, #interations until convergence, number of sites, number of samples at each site
                handle.write(cv.collapse_array_to_string(meta, str(i)))


def make_test_intervals(n):
    # hardcoded for now
    # more extreme cases maybe later
    unequal_splits = [[0.1, 0.9], [0.3, 0.7], [0.5, 0.5], [0.2, 0.2, 0.2, 0.2, 0.2], [0.1, 0.1, 0.2, 0.2, 0.4],[0.1, 0.1, 0.1, 0.1, 0.6], [0.2375, 0.2375, 0.2375, 0.2375, 0.05]]

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


if __name__ == "__main__":
    print('run split script')

    # parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    # parser.add_argument('-v', metavar='explained_var', type=float, help='explained variance')
    # parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    # parser.add_argument('-m', metavar='mult_dims_ret', type=str, help='comma separated list of intermediate dimensions', default = 1)
    # parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default = 100)
    # args = parser.parse_args()
    #
    # inputfile = args.f
    # outfile = args.o
    # exp_var = args.v
    # mult_dims_ret = args.m
    # sep = args.s
    # dims = args.d

    inputfile = '/home/anne/Downloads/mnist/flat.csv'
    outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/mnist/'
    exp_var = 0.5
    sep = ','
    mult_dims_ret = '0.75,5,10'
    dims = 10

    mult_dims_ret = cv.parse_array(mult_dims_ret)
    summaryfile = cv.make_eigenvector_path(outfile, path.basename(path.dirname(inputfile)) + '/' + str(exp_var))
    run_and_compare_unequal(inputfile, summaryfile, dims=dims, scale_unit=False, sep=sep, reported_angles=20,exp_var=exp_var, rownames=None, header=0)
