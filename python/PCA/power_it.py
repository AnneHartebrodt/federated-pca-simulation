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
from import_export.import_data import CustomDataImporter
import os. path as path
import convenience as cv


import scipy.sparse.linalg as lsa
import numpy as np
import scipy.linalg as la



class DistributedPowerIteration():

    def __init__(self, logger = None, file='/tmp/'):
        if logger is None:
            self.logger = logging.getLogger('default')
        else:
            self.logger = logger
        self.file = file
        self.cd = CustomDataImporter()


    def noise_variance_pow_it(self, epsilon, p, L, delta):
        '''

        Args:
            epsilon:
            p: iteration rank
            L: number of iterations
            delta:

        Returns: noise variance to be added during each iteration of the noisy power
                    method

        '''
        v = np.sqrt(4 * p * L * np.log(1 / delta)) * (1 / epsilon)
        return v

    def L_nr_it(self,eig_k, eig_q1, d):
        '''

        Args:
            eig_k: eigenvalue of target rank
            eig_q1: eigenvalue of intermediate iteration rank
            d: number of dimensions

        Returns: The required numbers of iteration to reach convergence

        '''
        L = np.ceil((eig_k / (eig_k - eig_q1))*np.log(d))
        return L

    def communication_overhead_M(self,eig_k, eig_q1, d, p, s):
        M = (eig_k / (eig_k - eig_q1)) * s * p * d * np.log(d)
        return M

    def coherence(self,eigenvectors, n):
        '''
        the matrix coherence is defined as follows:
        A=UEU_t (Eigenvalue decomposition)
        A elem mxm
        coh(V) m * ||U||inf
        Args:
            eigenvectors: Eigenvectors from the singular value decomposition as
                            np array
            n: Number of samples

        Returns: The coherence of the matrix

        '''

        m = np.max(eigenvectors)
        coher = n * m
        return coher

    def epsilon_bound(self,coher, d, L, eig_k, eig_q1, epsilon, delta, p):
        '''

        Args:
            coherence: Upper bound on matrix coherence (must be an estimate)
            d: number of dimensions
            L: number of iterations
            eig_k: eigenvalue of target rank
            eig_q1: eigenvalue of intermediate iteration rank
            epsilon:
            delta:
            p: number of dimensions of the data

        Returns: an upper bound for the overall privacy loss after the required
                number of iterations to reach convergence.

        '''
        v = self.noise_variance_pow_it(epsilon, p, L, delta)
        print(v)
        bound = v * np.sqrt(coher * np.log(d) * np.log(L)) / (eig_k - eig_q1)
        return bound

    def generate_random_gaussian(self,n, m, sigma):
        draws = n * m
        noise = sc.random.normal(0, sigma, draws)
        print('finished sampling')
        # make a matrix out of the noise
        noise.shape = (n, m)
        # transform matrix into s
        return noise

    def variance_explained(self, eigenvalues, perc = 0.5):
        total_variance = sum(eigenvalues)
        percentages = eigenvalues/total_variance
        p = 0
        sum_perc = 0
        while sum_perc<perc and p<len(eigenvalues):
            sum_perc = sum_perc+percentages[p]
            p = p+1
        return p

    def local_step(self, Xi, data):
        Yi = np.dot(data, Xi)
        return Yi

    def pooling_step(self, Yis, current):
        converged = False
        Yi = np.concatenate(Yis, axis=0)
        Yi, R = la.qr(Yi)

        if np.abs(np.dot(np.transpose(Yi), current)-Yi.shape[1])<0.01:
            converged=True
        return Yi, converged


    def power_method(self,data, sigma, L, p, noise=False):
        noise_norms = []
        # U,T,UT = lsa.svds(data, 1000)
        X_0 = self.generate_random_gaussian(data.shape[0], p, sigma)
        X_0, R = la.qr(X_0)
        # order eigenvectors from largest to smallest, to achive
        # ordered eigenvectors
        ord = np.argsort(np.diagonal(R))
        X_0 = np.flip(X_0[:,ord], axis=1)

        for i in range(L):
            # U1, E, UT1 = lsa.svds(X_0,1000)
            if noise:
                G = self.generate_random_gaussian(data.shape[0], p , np.max(X_0) * sigma)
                noise_norms.append(la.norm(G.flatten()))
                X_0, R = la.qr(np.dot(data, X_0[:,0:p]) + G)
            else:
                X_0, R =  la.qr(np.dot(data, X_0[:, 0:p]))
        eigenvals = self.eigenvalues(X_0, data)
        proj_global = sc.dot(data, X_0)

        return (proj_global, X_0[:, 0:p], eigenvals, noise_norms)

    # dimension n samples d dimensions
    # fix target rank k, intermediate rank 1 and iteration rank q
    # k>=q and 2q<=p<d

    def eigenvalue(self,A, v):
        Av = A.dot(v)
        return v.dot(Av)

    def eigenvalues(self, eigenvectors, cov):
        eigenvals = []
        for v in range(eigenvectors.shape[1]):
            eigenvals.append(self.eigenvalue(cov, eigenvectors[:, v]))
        return eigenvals

    def e_upper(self,eig_k, eig_q, r, d):
        e_upper = (eig_q / eig_k) * min(1 / np.log(eig_k / eig_q), 1 / np.log(r * d))
        return e_upper


    def theorem2_2(self, noise, eig_k, eig_q1, eigenvectors, p, q, r):
        '''
        These are the constraints that need to be fulfilled for theorem 2.2
        at every iteration
        Args:
            noise: The matrix of gaussian noise
            eig_k: The k'th eignevector
            eig_q1: The q+1th eigenvector
            eigenvectors: The matrix of the top q eigenvectors
            p: the iteration rank
            q: the intermediate rank
            r: some fixed constant r

        Returns: true, if the conditions are fulfilled, noise norm and

        '''
        bound = self.e_upper(eig_k, eig_q1, r, d)
        noise_norm = la.norm(noise.flatten())
        UqG_norm = la.norm(np.dot(eigenvectors, noise))
        o = bound * (eig_k - eig_q1)
        o2 = o*((np.sqrt(p)-np.sqrt(q- 1))/(r*np.sqrt(d)))
        return(5*noise_norm<=o and UqG_norm <=o2, noise_norm, UqG_norm)

    def assumed_noise(self, eig_k, eig_q1, e):
        return e*(eig_k-eig_q1)


    def determine_parameters(self, cov, epsilon, delta, n, var_exp=0.7, k=10, q1=20, p=40, r=1):
        r  = 1 # r is 'some parameter'
        d = cov.shape[1]
        U, E, UT = lsa.svds(cov, n)
        E, indz = cv.extract_eigenvals(E)
        U = np.flip(np.delete(U, indz, axis=1), axis=1)
        k, idx = self.variance_explained(E, var_exp)
        q1 = 2*k
        p = 4*k

        coher = self.coherence(U, d)
        L = self.L_nr_it(E[k], E[q1], n)
        noise_variance = self.noise_variance_pow_it(epsilon, p, L, delta)

        # The question is: are those two epsilons supposed to be the same
        # I don't think so. They just have the same name.
        # this is the epsilon for the noise
        eu = self.e_upper(E[k], E[q1 - 1], r, d)
        assumed_noise = self.assumed_noise(E[k], E[q1], eu)

        # this is the utility parameter
        bound = self.epsilon_bound(coher, d, L, E[k], E[q1], epsilon, delta, p)
        params = {'coher': coher, 'L':L, 'sigma':noise_variance, 'e_upper': eu, 'assumed_noise':assumed_noise,
                  'bound':bound, 'k':k, 'p':p, 'q1':q1}
        return params


if __name__ == '__main__':
    print('distributed power iteration library file')

