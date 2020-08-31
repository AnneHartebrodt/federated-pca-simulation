import numpy as np
import pandas as pd
import scipy as sc
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import shared_functions as s
import convenience as cv


def perform_SVD(data, t1=10):
    """
    Performs a singular value decomposition local data
    :param cov: A covariance matrix
    :param r: The number of top principal components to be considered
    :return: U_r*S_r (The product of the matrices taking the top r colums/rows)
    """

    U, S, UT, nd = svd_sub(data, t1)
    # In case we want to use more values in the approximation
    nd = min(nd, t1)
    R = np.zeros((nd, nd))
    np.fill_diagonal(R, S[0:nd])
    U_r = UT[0:nd, :]
    P = np.dot(np.sqrt(R), U_r)
    return P

def svd_sub(cov, ndims):
    # the sparse matrix version of svd has better memory requirements, while being a little
    # slower
    # covariance matrix is positive semi definite so SVD= Eigenvalue decomposition
    # print(nd)
    nd = min(cov.shape[1] - 1, ndims)
    V, S, W = sc.sparse.linalg.svds(cov, nd)
    # Sparse returns eigenvalues in ascending order
    S, indx = s.extract_eigenvals(S)
    W = np.flip(np.delete(W, indx, 0), axis=0)
    nd = min(nd, len(S))
    return V, S, W, nd




def aggregate_partial_SVDs(svds, t2=100):
    """
    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:
    """
    svds = np.concatenate(svds, axis=0)
    ndim = min(t2, svds.shape[0] - 1)
    U, S, UT, nd = svd_sub(svds, ndim)
    # nd = self.variance_explained(S, var_explained)
    UT = np.transpose(UT)
    return UT[:, 0:ndim], S[0:ndim]


def standalone_pca(data, ndims=100):
    """
    This function performs a standard principal component analysis
    :param data: input numpy data nxd with n the number of samples and d the number of variable
    :param ndims: Number of dimensions to project onto
    :return: the projected data, The first ndim eigenvectors and eigenvalues
    """
    cov = s.compute_cov(data)
    V, S, W, nd = svd_sub(cov, ndims)
    W = np.transpose(W)  # eigenvectors
    # create projection matrix by multiplying by the first nd eigenvectors
    proj_global = np.dot(data, W[:, 0:nd])
    return (proj_global, W[:, 0:nd], S[0:nd])


