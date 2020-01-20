import numpy as np
import pandas as pd
import scipy as sc
import random as r
import math as math
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import random as rnd
import os as os
import logging
import copy
import scipy.spatial.distance as d
import argparse as ap



class Distributed_DP_PCA():

    def __init__(self, logger = None, file='/tmp/'):
        if logger is None:
            self.logger = logging.getLogger('default')
        else:
            self.logger = logger
        self.file = file

    def compute_noisy_cov(self,original, epsilon0, delta0, nrSamples=0, nrSites=0, noise=True):
        """
        Compute a covariance matrix of a matrix and add epsilon, delta scaled
        iid gaussian noise to each entry of the matrix.

        Note: As currently implemented this is prone to attack due to the
        imperfect gaussian distribution (some values are more likely than others
        due to the underlying architechture) (Mironov et al.)

        :param original: The matrix for which the covariance is computed
        :param epsilon: epsilon as in (e,d)- differential privacy
        :param delta: delta as in (e,d)-differential privacy
        :param noise: True if noise should be added
        :return: the noisy covariance matrix
        """

        # compute covariance matrix
        # n = number of samples
        n = original.shape[0]
        # np.dot is matrix product for 2 dimensional arrays
        cov = (1 / (n-1)) * sc.dot(original.transpose(),original)
        print(cov.shape)
        pd.DataFrame(cov).to_csv(path_or_buf= self.file+'cov.tsv', sep='\t', header=None, index=False)
        #print(cov)
        if noise: # else return covariance matrix as is
            #print(psutil.virtual_memory())
            # variance of the added noise
            # in the case of unequal number of samples at each site the noise needs to
            # be scaled differently in order to achieve globally the same variance for the noise
            # this however comes with local DP with is weaker than epsilon 0.
            #epsiloni = (nrSamples/(nrSites*n))*epsilon0
            sigma = (1 / (n * epsilon0)) * sc.sqrt(2 * np.log(1.25 / delta0))
            # noise drawn from a 0 mean sigma variance normal distribution
            draws = cov.shape[1] * cov.shape[0]
            noise = sc.random.normal(0, sigma, draws)
            print('finished sampling')
            # make a matrix out of the noise
            noise.shape = cov.shape
            # transform matrix into symetric matrix
            # This approach is to memory inefficient, because an indicator matrix has to be held in memory
            # Looping is prohibitive for runtime reasons
            # A symetric matrix has to be added. Here first the upper, then the lower part
            # will be added
            # upper triangular matrix
            noise = sc.triu(noise)
            cov = cov + noise
            # make lower triangular matrix without diagonal
            noise = noise.transpose()
            noise = sc.tril(noise, -1)
            # add noise matrix and original matrix element-wise
            cov = cov + noise
            print(cov.shape)
        pd.DataFrame(cov).to_csv(path_or_buf=self.file+'noisy.tsv', sep='\t', header=None, index=False)

        return cov

    def perform_SVD(self,noisy_cov, ndims):
        """
        Performs a singular value decomposition of noisy_cov of the form
        A=USV.

        :param noisy_cov: A (noisy) covariance matrix
        :param r: The number of top principal components to be considered
        :return: U_r*S_r (The product of the matrices taking the top r colums/rows)
        """
        nd = min(noisy_cov.shape[1]-1, ndims)
        if(noisy_cov.shape[1]>10):
            print('Large covariance matrix: Using sparse PCA for performance')
            nd = min(noisy_cov.shape[1] - 1, ndims)
            U, S, UT = lsa.svds(noisy_cov, nd)
            # For some stupid reason sparse svd is returned in increasing order
            S = np.flip(S)
            UT = np.flip(UT, axis=0)
        else:
            nd = min(noisy_cov.shape[1], ndims)
            U, S, UT = la.svd(noisy_cov, lapack_driver='gesvd')

        R = np.zeros((nd, nd))
        np.fill_diagonal(R, S[0:nd])
        U_r = UT[0:nd,:]
        P = sc.dot(np.sqrt(R), U_r)
        return P

    def aggregate_partial_SVDs(self, svd_list, ndims=4, k=None):
        """
        This function aggregates the local proxy covariances by averaging them

        Function assumes equally shaped covariances matrices.
        :param svd_list: List of local P matrices
        :return:
        """
        #svd_list=self.normalize_eigenspaces(svd_list)
        # Average covariance matrices
        self.logger.info('Aggregating partial SVDs...')
        s = len(svd_list)

        Ac = np.dot(svd_list[0].transpose(), svd_list[0])
        for svd in range(1, len(svd_list)):
            Ac = Ac +np.dot(svd_list[svd].transpose(), svd_list[svd])
       # Ac = np.concatenate(svd_list)
        Ac = 1/s* Ac

        if (Ac.shape[1] > 10):
            print('Large covariance matrix: Using sparse PCA for performance')
            nd = min(Ac.shape[1] - 1, ndims)
            U, S, UT = lsa.svds(Ac, nd)
            # For some stupid reason sparse svd is returned in increasing order
            S = np.flip(S)
            UT = np.flip(UT, axis=0)
        else:
            U, S, UT = la.svd(Ac, lapack_driver='gesvd')

        UT = np.transpose(UT)
        UT = self.normalize_eigenvectors(UT)
        self.logger.info('...done')
        return UT[:, 0:ndims],S[0:ndims]

    def local_PCA(self, original, epsilon0, delta0, noise=True, ndims = 10):
        noisy_cov = self.compute_noisy_cov(original, epsilon0, delta0, noise=noise)
        PC = self.perform_SVD(noisy_cov, ndims)
        return PC

    def normalize_eigenvectors(self, V):
        """
        This function makes eigenvectors comparable, by assuring that the first element is
        positive and multipliing the vector by -1 otherwise.
        :param V: Eigenvector matrix with eigenvectors as column vectors
        :return: 'normalised' eigenvectors
        """
        for v in range(V.shape[1]):
            if V[0, v] <0:
                V[:, v]= V[:, v]*-1
        return V

    def normalize_eigenspaces(self, svd_list):
        for i in range(1, len(svd_list)):
            for column in range(svd_list[i].shape[1]):
                dp = sc.dot(svd_list[0][:,column], svd_list[i][:,column])
                if dp<0:
                    svd_list[i][:, column]= svd_list[i][:, column]*-1
                    print('norm')
        return svd_list



    def projection(self, scaled, sim, ndims, filename=None):
        projection = sc.dot(scaled, sim[:, 0:ndims])
        projection = pd.DataFrame(projection)
        if filename is not None:
            projection.to_csv(filename, sep='\t', header=None, index=False)
        return projection

    def calculate_euclidean_distance(self, V1, V2):
        res = []
        for line in range(min(V1.shape[1], V2.shape[1])):
            res.append(d.euclidean(V1[:, line], V2[:, line]))
        return (res)

    def standalone_pca(self, data, ndims):
        """
        This function performs a standard principal component analysis via eigendecomposition
        of the covariance matrix
        :param data: input numpy data nxd with n the number of samples and d the number of variable
        :param ndims: Number of dimensions to project onto
        :param center: if true, data is centered
        :param scale: if true, data is scaled to equal variance
        :param scale01: if true, data is scaled between 0 and 1
        :return: the projected data, The first ndim eigenvectors and eigenvalues
        """
        n = data.shape[0]  # get the number of rows
        cov = (1 / (n - 1)) * sc.dot(data.transpose(), data)  # use unbiased estimator
        # print(cov)
        # the sparse matrix version of svd has better memory requirements, while being a little
        # slower
        # covariance matrix is positive semi definite so SVD= Eigenvalue decomposition
        nd = min(data.shape[1], ndims)
        # print(nd)
        if (nd > 10):
            print('Using sparse PCA for decreased memory consumption')
            nd = min(data.shape[1] - 1, ndims)
            V_global, S, W = sc.sparse.linalg.svds(cov, nd)
            # For some stupid reason sparse svd is returned in increasing order
            S = np.flip(S)
            W = np.flip(W, axis=0)
        else:
            print('Using canonical PCA')
            nd = min(data.shape[1], ndims)
            V_global, S, W = sc.linalg.svd(cov)


        W = np.transpose(W)  # eigenvectors
        #W = self.normalize_eigenvectors(W)
        # this changes nothing

        # create projection matrix by multiplying by the first nd eigenvectors
        proj_global = sc.dot(data, W[:, 0:nd])
        return (proj_global, W[:, 0:nd], S[0:nd])



if __name__=="__main__":
    print('run')