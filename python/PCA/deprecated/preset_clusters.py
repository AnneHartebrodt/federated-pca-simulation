import argparse as ap
import os.path as path

import numpy as np
import pandas as pd

import proxy_covariance_runner as runner
import comparison as co
import convenience as cv
import proxy_covariance as dpca
import outlier_removal as orm
from import_export.import_data import CustomDataImporter


def run_and_compare_cluster_split(datafile, clusterfile, outfile=None, dims=1000, header=0, rownames=0, center=True,
                                  scale_var=True, scale01=False, scale_unit=False, transpose=False, sep='\t',
                                  reported_angles=20, exp_var=0.5):
    # ID is the folder name
    study_id = path.basename(path.dirname(datafile))
    print('Standalone PCA before outlier removal')
    importer = CustomDataImporter()
    data, varnames, sampleids = importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile,
                                                     transpose=False, sep=sep)
    n = data.shape[0]
    dims = min(dims, n)

    pca, W1, E1 = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames, center=center,
                                        scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                        transpose=transpose, sep=sep, log=True, exp_var=exp_var)

    print('Logging outliers')
    outliers = orm.outlier_removal_mad(pca, 6, 3)
    print(outliers)
    with open(outfile + '/removed_outliers.tsv', 'a+') as handle:
        handle.write(cv.collapse_array_to_string(outliers, study_id))

    print('Standalone PCA after outlier removal')
    pca, W1, E1 = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames, center=center,
                                        scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                        transpose=transpose, sep=sep, drop_samples=outliers, log=True, exp_var=exp_var)

    W, X = cluster_split(data, clusterfile=clusterfile, ndims=dims, exp_var=exp_var)
    id = study_id + '\t' + clusterfile + '\t' + str(exp_var)

    angles = co.compute_angles(W1, W, reported_angles=reported_angles)
    with open(outfile + '/angles_cluster_splits.tsv', 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=id))
    with open(outfile + '/eigenvalues.tsv', 'a+') as handle:
        handle.write(cv.collapse_array_to_string(X[0:reported_angles], study_id=id))


def cluster_split(data, clusterfile, scale_variance=True, center=True, scale01=False, scale_unit=False,
                  ndims=1000, header=0, rownames=0, scale_var=True, transpose=False, sep='\t', exp_var=0.5):
    """
    This function simulates a multisite PCA with each site having
    varying number of samples.
    :param data:
    :param sites: number of sites to split the data into
    :return:
    """

    Ac = []
    vexx = []
    s = 0
    start = 0
    importer = CustomDataImporter()
    clusters = pd.read_csv(filepath_or_buffer=clusterfile, header=header, sep=sep)

    for i in range(max(clusters.iloc[:, 1])):

        index = clusters.iloc[:, 1] == (i + 1)
        index = clusters[index].index
        # slice matrix
        data_sub = data[index, :]

        print('Local PCA for outlier identification')
        pca, W1, E1 = runner.run_standalone(data_sub, outfile, dims=ndims, header=header,
                                            rownames=rownames, center=center, scale_var=scale_var,
                                            scale01=scale01, scale_unit=scale_unit, transpose=transpose, sep=sep,
                                            drop_samples=[], log=True,
                                            exp_var=exp_var)

        print('Outlier removal')
        outliers = orm.outlier_removal_mad(pca, 6, 3)
        if len(outliers) != 0:
            data_sub = np.delete(data_sub, outliers, 0)

        print('Drop empty variable, scale local data, compute covariance')
        data_sub, vn = importer.drop0Columns(data_sub, None, drop=False, noise=True)
        data_sub = importer.log_transform(data_sub)
        data_sub = importer.scale_data(data_sub, center=center, scale_var=scale_variance, scale01=scale01,
                                       scale_unit=scale_unit)
        print('calculating cov')
        noisy_cov = dpca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, noise=False)
        print('finished')
        A, vex = dpca.perform_SVD(noisy_cov, ndims=ndims, mult_dims_returned=ndims, var_exp=exp_var)
        Ac.append(A)
        vexx.append(vex)

    W, X = dpca.aggregate_partial_SVDs(Ac, ndim=ndims)
    return W, X


if __name__ == "__main__":
    print('run split script')
    #
    # parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    # parser.add_argument('-v', metavar='explained_var', type=float, help='explained variance')
    # parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    # parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default=100)
    # parser.add_argument('-c', metavar='clusters', type=str, help='field delimiter', default=100)
    # args = parser.parse_args()
    #
    # inputfile = args.f
    # outfile = args.o
    # exp_var = args.v
    # sep = args.s
    # dims = args.d
    # clusterfile = args.c

    inputfile ='/home/anne/Documents/featurecloud/data/tcga/data_clean/CPTAC-2/coding_trunc.tsv'
    outfile = '/home/anne/Documents/featurecloud/results/test/target/'
    exp_var = 0.5
    sep = '\t'
    dims = 100
    clusterfile = '/home/anne/Documents/featurecloud/results/pca_plots/cluster/CPTAC-2_10_clusters.tsv'

    run_and_compare_cluster_split(inputfile, clusterfile, outfile=outfile, dims=dims, header=0, rownames=0, center=True,
                                  scale_var=True, scale01=False, scale_unit=False, transpose=False, sep='\t',
                                  reported_angles=20, exp_var=exp_var)
