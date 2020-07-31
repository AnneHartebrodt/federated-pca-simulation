import os.path as path

import math as math
import numpy as np
import scipy.linalg as la
import scipy.linalg as lsa

import fast_power as powerit
import import_export.easy_import as easy
import import_export.import_data as imp
import guo_vertical as gv


def run_standalone(data, outfile=None, p=10, header=None, rownames=None, center=False, scale_var=False, scale01=False,
                   scale_unit=False, transpose=False, sep='\t', drop_samples=[], log=False, nr_iter=10, qr=True):
    """
    This function performs a regular principal component analysis and saves the result to files containing
    the projection the
    :param datafile: Unscaled datafile
    :param ddpca:
    :param outfile: path and name of the output file without extension
    :param dims: Number of dimensions to return (#eigenvectors and corresponding eigenvalues)
    :param seed: random seed
    :param nr_samples: #variables to select if not all the data columns are to be used for the pca
    :param header: row number which contains the header/ number of header rows
    :param rownames: column number which contains the rownames/sample ids
    :return: projection, eigenvectors and eigenvalues
    """
    # if data is a string, data needs to be read first, otherwise it is
    # assumed to be scaled and ready to use
    if isinstance(data, str):
        data = easy.easy_import(data, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,
                                scale01=scale01, scale_unit=scale_unit, drop_samples=drop_samples, log=log,
                                outfile=outfile, transpose=transpose)
    # run pca
    if qr:
        UT, ST, VT = powerit.power_method_via_qr(data=data, p=p, maxit=nr_iter)
    else:
        UT, ST, VT = powerit.power_method(data=data, p=p, maxit=nr_iter)
    # return the projected datapoints, the eigenvectors and the eigenvalues
    return UT, ST, VT


def simulate_distributed(local, p=10, maxit=10):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """

    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors

    locals = []
    for l in local:
        locals.append(powerit.local_step(data=l, maxit=maxit, p=p))

    U = powerit.pooling_step(locals, p=p - 1)
    loc_ev = []
    for l in local:
        loc_ev.append(powerit.local_step_2(data=l, U=U, p=p - 2))

    return loc_ev, U


def simulate_distributed_2(local, p=8, maxit=10):
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
    G_is = []
    for i in range(len(local)):
        # initialise matrices randomly
        G_is.append(powerit.generate_random_gaussian(d, p))
    H_is = []

    maxit = max(maxit, 1)
    for i in range(maxit):
        for l, G_i in zip(local, G_is):
            # calculate H matrices
            H_is.append(powerit.local_1(l, G_i))
        # combine Hi matrices to one
        H_i = powerit.pooling_step(H_is)
        H_is = []
        # delete G_is and replace by new guess
        G_is = []
        for l in local:
            G_is.append(powerit.local_2(H_i, l))

    loc_ev = []
    for l in local:
        loc_ev.append(powerit.local_step_2(data=l, U=H_i))
    return loc_ev, H_i


def simulate_distributed_multi_local(local, p=8, weights=None, nr_local_rounds=5, tolerance=0.000001):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[0]
    X = powerit.generate_random_gaussian(d, p, sigma=1)
    X, R = lsa.qr(X, mode='economic')
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    while not converged:
        count = count + 1
        locals = []
        for l in local:
            # run a few iterations locally before sending the result to the aggrgator
            for i in range(nr_local_rounds):
                locals.append(powerit.local_step(X, l))
        X, E, converged = powerit.pooling_step(locals, X, weights=weights)
    ord = np.argsort(E)
    X = np.flip(X[:, ord], axis=1)
    E = np.flip(np.sort(E))
    return X, E, count


def chunk_data(data, interval_end):
    data = data.T
    np.random.shuffle(data)
    localdata = []
    interval_end = [np.floor(data.shape[0] * x) for x in interval_end]
    start = 0

    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[start:end, :]
        # calculate covariance matrix
        start = int(interval_end[i])
        localdata.append(data_sub.T)
    return localdata


def angleo(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree
    Args:
        v1: first vector
        v2: second vector

    Returns: angle in degree or NaN

    """
    dp = np.dot(v1, v2)
    norm_x = la.norm(v1)
    norm_y = la.norm(v2)
    co = np.clip(dp / (norm_x * norm_y), -1, 1)
    theta = np.arccos(co)
    a = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(a):
        return a
    # return the canonical angle
    return a


def eigenvector_assembler(eigenvectors):
    E = []
    for e in eigenvectors:
        a = angleo(e.T, np.random.random(len(e)))
        if a < 90:
            print(str(a) + 'notflipped')
            E = np.concatenate([E, e])
        else:
            E = np.concatenate([E, e * -1])
            print(str(a) + 'flipped')
    E = E / np.linalg.norm(E)
    return E


def eigenvector_assembler_2_dumb(eigenvectors):
    E = np.concatenate([eigenvectors[0], eigenvectors[1]])
    E2 = np.concatenate([eigenvectors[0], -1 * eigenvectors[1]])
    return E, E2


def mev(u, truth):
    k = min(truth.shape[1], u.shape[1])  # number of eigenvectors in subspace
    print(k)
    m = np.dot(u.T, truth)
    sum = 0
    for i in range(min(truth.shape[1], u.shape[1])):
        sum = sum + np.linalg.norm(m[:, i], 2)
    mev = sum / k
    return mev


def make_result(eigenvectors, v):
    res = []
    for i in range(v.shape[1]):
        res.append(co.angle(eigenvectors[:, i], v[:, i]))
    for i in range(v.shape[1]):
        res.append(mev(eigenvectors[:, 0:i + 1], v[:, 0:i + 1]))
    return res


def subspace_reconstruction_error(data, eigenvectors):
    res = []
    for i in range(v.shape[1]):
        proj = np.dot(data, eigenvectors[:, 0:i])
        rec = np.dot(proj, eigenvectors[:, 0:i].T)
        res.append(np.linalg.norm(data - rec) / (data.shape[1] * data.shape[0]))
    return res


if __name__ == '__main__':

    import scipy.sparse.linalg as lsa
    import matplotlib.pyplot as plt
    import comparison as co

    file = '/home/anne/Documents/data/radiation/data_new.tsv'
    outfile = '/home/anne/Documents/featurecloud/gwas/chr10'
    header = 0
    rownames = None
    center = False
    scale = False
    scale_var = False
    scale01 = False
    scale_unit = False
    p = 23
    transpose = True
    #
    import pandas as pd

    # import the eigenvectors produced by plink2 as a comparison
    eigenvector_plink = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/1000g10/plink2.eigenvec', header=0)
    eigenvector_plink = np.array(eigenvector_plink)
    eigenvector_plink = eigenvector_plink[:, 2:eigenvector_plink.shape[1]]

    eigenvalues_plink = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/1000g10/plink2.eigenval',
                                      header=None)
    eigenvalues_plink = np.asarray(eigenvalues_plink)
    # import and center data, because this is apparently done automaticall in fastPCA
    data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                            scale01=scale01, scale_unit=scale_unit,
                            outfile=outfile, transpose=transpose, sep='\t')
    data = imp.CustomDataImporter().scale_center_data(data, center=True)

    # regular PCA as implemented in scipy
    # as we use the sparse version of scipy we need to flip the order of the eigenvectors

    final_dim = 10
    u, s, v = lsa.svds(data, k=final_dim)
    u = np.flip(u, axis=1)
    v = np.flip(v.T, axis=1)
    s = np.flip(s)

    nr_repeats = 10
    # make a results dataframe containing
    # for both qr and
    # - subspace reconstruction error (1)
    # - mev (final_dim)
    # - angles between eigenvectors (final_dim)
    # - eigenvalues (final_dim)

    eigenvalues_qr = []
    eigenvalues_svd = []

    subspace_rec_err_qr = []
    subspace_rec_err_svd = []

    results_qr = []
    results_svd = []

    results_qr_vs_plink = []
    results_svd_vs_plink = []
    # repeat the standalone PCA to be sure the result is not influenced by random init
    for i in range(nr_repeats):
        UT_qr, ST_qr, VT_qr = run_standalone(data, outfile=outfile, p=p, header=header, rownames=rownames,
                                             center=center, scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                             transpose=transpose, nr_iter=20, qr=True)

        UT_svd, ST_svd, VT_svd = run_standalone(data, outfile=outfile, p=p, header=header, rownames=rownames,
                                                center=center, scale_var=scale_var, scale01=scale01,
                                                scale_unit=scale_unit,
                                                transpose=transpose, nr_iter=20, qr=False)

        results_qr.append(make_result(VT_qr, v))
        results_svd.append(make_result(VT_svd, v))

        # results_qr_vs_plink.append(make_result(VT_qr, eigenvector_plink))
        # results_svd_vs_plink.append(make_result(VT_svd,eigenvector_plink))

        eigenvalues_qr.append(np.sqrt(ST_qr))
        eigenvalues_svd.append(np.sqrt(ST_svd))

        subspace_rec_err_qr.append(subspace_reconstruction_error(data, VT_qr))
        subspace_rec_err_svd.append(subspace_reconstruction_error(data, VT_svd))

    chunks = [0.5, 1.0]
    result_fed = []
    result_fed_2 = []
    result_fed_vs_plink = []
    eigenvalues_fed = []
    subspace_rec_fed = []

    for dim in range(nr_repeats):
        loca = chunk_data(data, chunks)
        loc_ev, H_i = simulate_distributed(loca, maxit=20, p=p)

        l = {}
        eigenvalues = []
        eigenvectors_assembled = []
        eigenvector_ass2 = []
        for i in range(len(loc_ev)):
            eigenvalues.append(loc_ev[i][1])
            for e in range(loc_ev[i][2].shape[1]):
                if str(e) in l:
                    l[str(e)].append(loc_ev[i][2][:, e])
                else:
                    l[str(e)] = []
                    l[str(e)].append(loc_ev[i][2][:, e])

        eigenvalues = np.asarray(eigenvalues)
        eigenvalues = np.mean(eigenvalues, axis=0)
        eigenvalues = np.sqrt(eigenvalues)
        eigenvalues_fed.append(eigenvalues)

        for e in l.keys():
            if int(e) < v.shape[1]:
                c, c2 = eigenvector_assembler_2_dumb(l[e])
                eigenvectors_assembled.append(c)
                eigenvector_ass2.append(c2)

        eigenvectors_assembled = np.asarray(eigenvectors_assembled).T
        eigenvector_ass2 = np.asarray(eigenvector_ass2).T
        result_fed.append(make_result(eigenvectors_assembled, v))
        result_fed_2.append(make_result(eigenvector_ass2, v))
        # result_fed_vs_plink.append(make_result(eigenvectors_assembled, eigenvector_plink))

        subspace_rec_fed.append(subspace_reconstruction_error(data, eigenvectors_assembled))

    pd.DataFrame(eigenvalues_qr).to_csv(path.join(outfile, 'eigenvalues_qr.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(eigenvalues_svd).to_csv(path.join(outfile, 'eigenvalues_svd.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(eigenvalues_fed).to_csv(path.join(outfile, 'eigenvalues_fed.tsv'), sep='\t', header=None, index=None)

    pd.DataFrame(subspace_rec_err_qr).to_csv(path.join(outfile, 'sub_rec_err_qr.tsv'), sep='\t', header=None,
                                             index=None)
    pd.DataFrame(subspace_rec_err_svd).to_csv(path.join(outfile, 'sub_rec_err_svd.tsv'), sep='\t', header=None,
                                              index=None)
    pd.DataFrame(subspace_rec_fed).to_csv(path.join(outfile, 'sub_rec_err_fed.tsv'), sep='\t', header=None, index=None)

    result_svd = np.asarray(results_svd)
    result_qr = np.asarray(results_qr)
    result_fed = np.asarray(result_fed)
    result_qr_vs_plink = np.asarray(results_qr_vs_plink)
    result_svd_vs_plink = np.asarray(results_svd_vs_plink)
    result_fed_vs_plink = np.asarray(result_fed_vs_plink)

    pd.DataFrame(result_fed_vs_plink).to_csv(path.join(outfile, 'fed_vs_plink.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(result_svd_vs_plink).to_csv(path.join(outfile, 'svd_vs_plink.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(result_qr_vs_plink).to_csv(path.join(outfile, 'qr_vs_plink.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(result_qr).to_csv(path.join(outfile, 'scipy_vs_qr.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(result_svd).to_csv(path.join(outfile, 'scipy_vs_svd.tsv'), sep='\t', header=None, index=None)
    pd.DataFrame(result_fed).to_csv(path.join(outfile, 'scipy_vs_fed.tsv'), sep='\t', header=None, index=None)

    proj = np.dot(data, VT_qr)
    plt.scatter(proj[:, 0], proj[:, 1])
    plt.show()

    plt.scatter(proj[:, 1], proj[:, 2])
    plt.show()

    plt.scatter(proj[:, 2], proj[:, 3])
    plt.show()

    result_fed_2 = np.asarray(result_fed_2)

    # print('power iteration runner')
    # parser = ap.ArgumentParser(description='Run distributed PCA simulation')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-d', metavar='dimensions', type=int, help='number of principal components to return')
    # parser.add_argument('-p', metavar='output directory', type=str, help='output directory for simulation study.')
    # parser.add_argument('-k', metavar='number_hospitals', type=int, help='Number of simulated hospitals', default=5)
    # parser.add_argument('-s', action='store_true',
    #                     help='If true the generated eigenspaces are saved (!a lot of large files!)', default=False)
    # parser.add_argument('-c', help='True of data has column headers',
    #                     default=False, action='store_true')
    # parser.add_argument('-i', metavar='sampleids', type=int,
    #                     help='Dataframe column which contains the sample ids, 0 index,', default=-1)
    # parser.add_argument('-v', action='store_true',
    #                     help='Scale variables to unit variance', default=True)
    # parser.add_argument('-u', action='store_true',
    #                     help='Scale samples to unit norm', default=True)
    # parser.add_argument('-z', action='store_true',
    #                     help='Scale variables between 0 and 1', default=False)
    # parser.add_argument('-t', action='store_true',
    #                     help='Center variables by substracting the mean', default=True)
    #
    # parser.add_argument('-A', action='store_true',
    #                     help='Run standalone simulation', default=False)
    # args = parser.parse_args()
    #
    #
    # def parse_array(value_str):
    #     values_str = value_str.split(',')
    #     values_int = []
    #     for v in values_str:
    #         values_int.append(float(v))
    #     return values_int
    #
    #
    # if args.c:
    #     header = 0
    # else:
    #     header = None
    #
    # if args.i == -1:
    #     rownames = None
    # else:
    #     rownames = args.i
    #
    # try:
    #     epsilons = parse_array(args.e)
    # except:
    #     print('Incompatible epsilon parameters')
    #     print('default to epsilon=0.1')
    #     epsilons = [0.1]
    #
    # try:
    #     deltas = parse_array(args.g)
    # except:
    #     print('Incompatible delta parameters')
    #     print('default to delta=0.001')
    #     deltas = [0.001]
    #
    # print('file: ' + str(args.f))
    # print('dimensions: ' + str(args.d))
    # print('output directory: ' + str(args.p))
    # print('number of hospitals: ' + str(args.k))
    # print('save eigenvalues: ' + str(args.s))
    # print('repeats: ' + str(args.r))
    # print('column headers: ' + str(args.c))
    # print('sample ids: ' + str(args.i))
    #
    # if args.A:
    #     run_standalone(args.f, outfile=args.p, p=args.d, header=header, rownames=rownames,
    #                    center=args.t, scale_var=args.v, scale01=args.z, scale_unit=args.u,
    #                    transpose=False)
