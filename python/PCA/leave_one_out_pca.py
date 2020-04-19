import argparse as ap
import itertools as it
import os as os
import os.path as path

import numpy as np
import pandas as pd

import proxy_covariance_runner as runner
import comparison as co
import convenience as cv
import outlier_removal as orm
from import_export.import_data import CustomDataImporter


def superset(set):
    la = []
    if len(set) < 5:
        ll = list(it.product([0, 1], repeat=len(set)))
    else:
        ll = []
        for i in range(20):
            ll.append(np.random.choice([0, 1], size=len(set)))
    for l in range(len(ll)):
        tmp = []
        for i in range(len(ll[l])):
            if ll[l][i] == 1:
                tmp.append(set[i])
        if len(tmp) > 1 and len(tmp) < len(set):
            la.append(tmp)
    for s in set:
        la.append([s])
    return la


def run_dropout(data_scaled, outdir, tempdir, nr_dropped, study_id, outliers=[], mode='outlier_free', ndims=100,
                header=None, rownames=None, center=True, scale_var=True, scale01=False, scale_unit=False,
                transpose=False, sep='\t'):
    '''
    This function runs a 'leave-k'out' simulation. An eigenvalue decomposition
    is calculated n times on a datamatrix, leaving out nr_dropped samples at each time.
    The indices of the dropped samples, the eigenvalues and the eigenvectors are saved
    (up to 20).
    Args:
        data_scaled: The scaled data matrix
        nr_dropped: The number of samples to drop at each pass

    Returns: None, the according results are saved in files

    '''

    n = data_scaled.shape[0]

    # Outlier free mode means, we can estimate the influence of a sample
    # under ideal conditions
    if mode == 'outlier_free':
        if len(outliers) != 0:
            data_scaled = np.delete(data_scaled, outliers, 0)

    data_scaled_copy = np.copy(data_scaled)
    # the number of draws is 20% of the number of samples,
    # but a minimum of 10 samples, unless the dataframe is smaller
    # in which case very sample is dropped once.

    draws = max(int(np.ceil(n / 5)), min(10, n))

    if mode == 'cheat':
        # gene
        droprows = superset(outliers)
        nr_dropped = 1
    else:
        droprows = np.random.choice(n, size=draws, replace=False)
    index = 0

    #
    for row in range(0, len(droprows), nr_dropped):
        if mode == 'cheat':
            data_scaled = np.delete(data_scaled, droprows[row], 0)
        else:
            # remove the samples rows= patients from the data
            e = min(row + nr_dropped, data_scaled.shape[0])
            data_scaled = np.delete(data_scaled, droprows[row:e], 0)
        # calculate covariance matrix and svd
        # print('Computing Singular Value Decomposition')
        # It should be sufficient to calculate at most n PCs, because the others will contain noise
        nd = min(min(n - nr_dropped, data_scaled.shape[0] - 1), ndims)
        pca, U, E = runner.run_standalone(data_scaled, outdir, dims=ndims, header=header, rownames=rownames,
                                          center=center, scale_var=scale_var, scale01=scale01,
                                          scale_unit=scale_unit, transpose=transpose, sep=sep, log=True)
        E = E[0:20]
        U = U[:, 0:20]

        # print the 20 first eigenvectors
        res_eigen = cv.collapse_array_to_string_nrdropped(E, study_id, nr_dropped)
        # save the dropped patients, to trace back big changes later
        res_dropped = cv.collapse_array_to_string_nrdropped(droprows[row:(row + nr_dropped)], study_id, nr_dropped)

        # write eigenvalues, dropped samples and eigenvectors to file
        with open(path.join(outfile, 'eigenvalues.tsv'), 'a+') as handle:
            handle.write(res_eigen)
        with open(path.join(outfile, 'dropped_samples.tsv'), 'a+') as handle:
            handle.write(res_dropped)
        # save the eigenvectors in a separate directory. They likely will be deleted after
        # calculating the angles due to missing diskspace
        pd.DataFrame(U).to_csv(path.join(tempdir, 'eigenvectors' + str(index) + '.tsv'), sep='\t', header=None, index=None)
        # increment the index
        index = index + 1
        # get the correct data back.
        data_scaled = np.copy(data_scaled_copy)


def run_standalone(data, outdir=None, dims=100, header=None, rownames=None, center=True, scale_var=True, scale01=False,
                   scale_unit=False, transpose=False, sep='\t', reported_angles=20, exp_var=1, study_id='ID1'):
    pca, W1, E1 = runner.run_standalone(data, outdir, dims=dims, header=header, rownames=rownames,
                                        center=center, scale_var=scale_var, scale01=scale01,
                                        scale_unit=scale_unit, transpose=transpose, sep=sep, log=True)

    res_eigen = cv.collapse_array_to_string_nrdropped(E1[0:reported_angles], study_id, '0')
    with open(path.join(outfile, 'eigenvalues_with_outliers.tsv'), 'a+') as handle:
        handle.write(res_eigen)
    # pd.DataFrame(W1[:,0:20]).to_csv(ev_path + '/eigenvectors_with_outlier.tsv', sep='\t', header=None, index=None)

    print('Logging outliers')
    outliers = orm.outlier_removal_mad(pca, 6, 3)
    print(outliers)
    with open(path.join(outfile, 'removed_outliers.tsv'), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(outliers, study_id))

    # delete outliers manually
    if len(outliers) != 0:
        data = np.delete(data, outliers, 0)

    print('Standalone PCA after outlier removal')
    pca, W1, E1 = runner.run_standalone(data, outdir, dims=dims, header=header, rownames=rownames,
                                        center=center, scale_var=scale_var, scale01=scale01,
                                        scale_unit=scale_unit, transpose=transpose, sep=sep, log=True)

    res_eigen = cv.collapse_array_to_string_nrdropped(E1[0:reported_angles], study_id, '0')

    # write eigenvalues, dropped samples and eigenvectors to file
    with open(path.join(outfile, 'eigenvalues_reference.tsv'), 'a+') as handle:
        handle.write(res_eigen)
    pd.DataFrame(W1[:, 0:reported_angles]).to_csv(path.join(outdir, 'eigenvectors_reference.tsv'), sep='\t', header=None,
                                                  index=None)
    return outliers


def run_study(datafile, outdir, tempdir, header=None, rownames=None, sep='\t', dims=100, mode='outlier_free',
              nr_dropped=1):
    # ID is the folder name
    study_id = path.basename(path.dirname(datafile))

    importer = CustomDataImporter()
    data, varnames, sampleids = importer.data_import(datafile, header=header, rownames=rownames, outfile=outdir,
                                                     transpose=False, sep=sep)
    dims = min(dims, data.shape[0])
    outliers = run_standalone(data, outdir=outdir, study_id=study_id)
    run_dropout(data, outdir, tempdir, nr_dropped=nr_dropped, study_id=study_id, outliers=outliers, mode=mode,
                ndims=dims)
    compute_angles(tempdir, outdir, nr_dropped)


def compute_angles(tempdir, outdir, nr_dropped):
    '''
    computes the angle in degrees for all the input files
    in the given directory vs the canoncial pca. Angles are between the matching vectors
    (columnwise)

    Args:
        eigenvector_path: folder containing files, which contain
        columnwise eigenvectors in decreaseing order according to the
        eigenvector
        other_folder: a different folder, where the outputs are saved.
        Make sure it is a different folder in case the eigenvector folder
        is deleted.

    Returns:

    '''

    # read reference
    reference = pd.read_csv(filepath_or_buffer=path.join(outdir, 'eigenvectors_reference.tsv'), sep='\t', header=None,
                            engine='python')
    angles = []
    i = 1
    for file in os.listdir(tempdir):
        # read one of the dropped out eigenvalue files
        current = pd.read_csv(filepath_or_buffer=path.join(tempdir,file), sep='\t', header=None, engine='python')
        # calculate angles against reference
        for k in range(max(current.shape[1], reference.shape[1])):
            ll = [i, 0, k, co.angle(current.iloc[:, k], reference.iloc[:, k])]
            angles.append(ll)
        i = i + 1
    pd.DataFrame(angles).to_csv(path.join(outdir, 'angles_dropout' + str(nr_dropped) + '.tsv'), sep='\t', header=None,
                                index=None)


if __name__ == "__main__":
    print('run eigengap script')

    parser = ap.ArgumentParser(description='Eigengap calculation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default=100)
    parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    parser.add_argument('-n', metavar='nr_dropped', type=int, help='number of rows dropped')
    parser.add_argument('-m', metavar='mode', type=str, help='cheat or outflier_free')
    args = parser.parse_args()
    #
    inputfile = args.f
    outfile = args.o
    dims = args.d
    sep = args.s
    nr_dropped = args.n
    mode = args.m


    dname = path.join(path.basename(path.dirname(inputfile)), mode + '_' + str(nr_dropped))
    summaryfile = cv.make_eigenvector_path(outfile, dname)
    ev_path = cv.make_eigenvector_path(summaryfile, 'eigenvectors')
    run_study(inputfile, summaryfile, ev_path, sep=sep, header=0, dims=dims, mode=mode, nr_dropped=nr_dropped)
    cv.delete(ev_path)
