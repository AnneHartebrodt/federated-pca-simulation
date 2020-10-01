import numpy as np
import scipy as sc
import scipy.sparse.linalg as lsa
import pandas


def svd_sub(cov, ndims):
    # the sparse matrix version of svd has better memory requirements, while being a little
    # slower
    # covariance matrix is positive semi definite so SVD= Eigenvalue decomposition
    # print(nd)
    nd = min(cov.shape[1] - 1, ndims)
    V, S, W = sc.sparse.linalg.svds(cov, nd)
    # Sparse returns eigenvalues in ascending order
    S, indx = extract_eigenvals(S)
    W = np.flip(np.delete(W, indx, 0), axis=0)
    nd = min(nd, len(S))
    return V, S, W, nd


def compute_cov(original):
    """
    Compute a covariance matrix of a matrix
    :return: the noisy covariance matrix
    """
    # compute covariance matrix
    # n = number of samples
    n = original.shape[0]
    # np.dot is matrix product for 2 dimensional arrays
    cov = (1 / (n - 1)) * np.dot(original.transpose(), original)
    return cov


def variance_explained(eigenvalues, perc=0.5):
    total_variance = sum(eigenvalues)
    percentages = eigenvalues / total_variance
    p = 0
    sum_perc = 0
    while sum_perc < perc and p < len(eigenvalues):
        sum_perc = sum_perc + percentages[p]
        p = p + 1
    return p


def generate_random_gaussian(n, m):
    draws = n * m
    noise = sc.random.normal(0, 1, draws)
    print('Generated random intial matrix: finished sampling')
    # make a matrix out of the noise
    noise.shape = (n, m)
    # transform matrix into s
    return noise


def extract_eigenvals(eigenvalues):
    """
    eigenvaluesigendecomposition from scipy.linalg.sparse returns eigenvalues ordered in
    increasing order, followed by eigenvalues which are 0.
    Eigenvalues are returned in decreasing order ommiting th 0s alltogether
    Args:
        eigenvalues: Eigenvalues from a sparse singular value decomposition.

    Returns: Eigenvalue vector in decreasing order, without 0s.

    """
    indz = np.where(eigenvalues == 0)
    eigenvalues = np.flip(eigenvalues)
    eigenvalues = eigenvalues[eigenvalues != 0]
    return eigenvalues, indz


def projection(scaled, sim, ndims):
    projection = sc.dot(scaled, sim[:, 0:ndims])
    return projection


def partition_data_horizontally(data, splits=2, equal=True):
    n = data.shape[0]
    interval_end = []
    if equal:
        for s in range(splits - 1):
            interval_end.append(int(np.floor((s + 1) * n / splits)))
        # last one should include all the elements!
        interval_end.append(n)

    np.random.shuffle(data)
    start = 0
    localdata = []
    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[start:end, :]
        # calculate covariance matrix
        start = int(interval_end[i])
        localdata.append(data_sub)
    return localdata


def partition_data_vertically(data, splits=2, equal=True, randomize=False):
    n = data.shape[1]
    interval_end = []
    if equal:
        for s in range(splits - 1):
            interval_end.append(int(np.floor((s + 1) * n / splits)))
        # last one should include all the elements!
        interval_end.append(n)
    print(n)
    print(n)
    if randomize:
        # generate a permutation of the indices and return the permuted data
        choice = np.random.choice(n, n, replace=False)
        data = data[:, choice]
    else:
        choice = range(n)
    start = 0
    localdata = []
    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[:, start:end]
        start = int(interval_end[i])
        localdata.append(data_sub)
    return localdata, choice


def eigenvalue(A, v):
    Av = A.dot(v)
    return v.dot(Av)


def eigenvalues(eigenvectors, cov):
    eigenvals = []
    for v in range(eigenvectors.shape[1]):
        eigenvals.append(eigenvalue(cov, eigenvectors[:, v]))
    return eigenvals
