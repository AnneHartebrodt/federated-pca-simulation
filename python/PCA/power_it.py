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
        L = np.ceil((eig_k * np.log(d)) / (eig_k - eig_q1))
        return L

    def communication_overhead_M(self,eig_k, eig_q1, d, p, s):
        M = (eig_k / (eig_k - eig_q1)) * s * p * d * np.log(d)
        return M

    def coherence(self,eigenvectors, n):
        '''
        the matrix coherence is defined as follows:
        A=UEU_t (Eigenvalue decomposition)
        A elem mxn
        coh(V) = max {m ||U||inf^2, n ||V||inf^2}
        Args:
            eigenvectors: Eigenvectors from the singular value decomposition as
                            np array
            n: Number of samples

        Returns: The coherence of the matrix

        '''

        m = np.max(eigenvectors)
        coher = n * m * m
        return coher

    def epsilon_bound(self,coherence, d, L, eig_k, eig_q1, epsilon, delta, p):
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
        bound = v * np.sqrt(coherence * np.log(d) * np.log(L)) / (eig_k - eig_q1)
        return bound

    def generate_random_gaussian(self,n, m, sigma):
        draws = n * m
        noise = sc.random.normal(0, sigma, draws)
        print('finished sampling')
        # make a matrix out of the noise
        noise.shape = (n, m)
        # transform matrix into s
        return noise

    def power_method(self,data, sigma, L, q, noise=False):
        noise_norms = []
        # U,T,UT = lsa.svds(data, 1000)
        X_0 = self.generate_random_gaussian(data.shape[0], q, sigma)
        X_0, R = la.qr(X_0)
        # order eigenvectors from largest to smallest, to achive
        # ordered eigenvectors
        ord = np.argsort(np.diagonal(R))
        X_0 = np.flip(X_0[:,ord], axis=1)

        for i in range(L):
            # U1, E, UT1 = lsa.svds(X_0,1000)
            if noise:
                G = self.generate_random_gaussian(data.shape[0], q , np.max(X_0) * sigma)
                noise_norms.append(la.norm(G.flatten()))
                X_0, R = la.qr(np.dot(data, X_0[:,0:q]) + G)
            else:
                X_0, R =  la.qr(np.dot(data, X_0[:, 0:q]))
        return (X_0[:, 0:q], noise_norms)

    # dimension n samples d dimensions
    # fix target rank k, intermediate rank 1 and iteration rank q
    # k>=q and 2q<=p<d

    def eigenvalue(self,A, v):
        Av = A.dot(v)
        return v.dot(Av)

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

    def extract_eigenvals(self, E):
        '''
        Eigendecomposition from scipy.linalg.sparse returns eigenvalues ordered in
        increasing order, followed by eigenvalues which are 0.
        Eigenvalues are returned in decreasing order ommiting th 0s alltogether
        Args:
            E: Eigenvalues from a sparse singular value decomposition.

        Returns: Eigenvalue vector in decreasing order, without 0s.

        '''
        E = np.flip(E)
        E = E[E != 0]
        return E


    def generate_result_str(self, inputfile, epsilon, delta, n, d, q, p, r, coherence, L, noise_variance, e, noise_vals):
        res = path.basename(path.dirname(inputfile)) + '\t' + str(epsilon)+ '\t' + str(delta) + '\t' + str(n)+'\t' + str(d)+'\t' + str(q)+'\t' + str(p)\
              +'\t'+str(r)+'\t'+ str(coherence)+ str(L)+'\t' + str(noise_variance)+ '\t'+ str(e)+'\t'
        header = 'study.id\tepsilon\tdelta\tnr.samples\tnr.dimension\tnr.itermediate\tit.rank\tr\tcoherence\tnr.iterations\tnoise.variance\terror.bound\t'
        i = 1
        for v in noise_vals:
            res = res +str(v) + '\t'
            header = 'n1.'+str(i)+'\t'
            i= i+1
        return (res, header)

    def write_summary(self, res, header, outfile):
        try:
            os.makedirs(path.dirname(outfile))
        except OSError:
            print(path.dirname(outfile) + ' could not be created')
        else:
            print(path.dirname(outfile) + ' was created')

        if not path.exists(outfile):
            with open(outfile, 'w') as handle:
                handle.write(header + '\n')

        with open(outfile, 'a+') as handle:
            handle.write(res + '\n')

if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Eigengap calculation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    args = parser.parse_args()

    inputfile = args.f
    outfile = args.o

    #inputfile = '/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/output_transposed.txt'
    #outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/summary.txt'


    cd = CustomDataImporter()
    DPIT = DistributedPowerIteration()

    # Scale data an calculate eigengap
    data, sample_id, var_names = cd.data_import(inputfile, header=0, rownames=0)
    data_scaled = cd.scale_data(data, center=True, scale_var=True, scale_unit=True)
    n = data.shape[0]

    cov = np.cov(data_scaled.T)
    U, E, UT = lsa.svds(cov, n)
    E = DPIT.extract_eigenvals(E)

    pd.DataFrame(E[0:15]).to_csv(path.dirname(inputfile)+'/eigenvalues.tsv', sep='\t')
    pd.DataFrame(np.flip(UT, axis=1)[0:15]).to_csv(path.dirname(inputfile)+'/eigenvectors.tsv', sep='\t')

    d = cov.shape[1]
    k = 5
    q1 = 15
    p = 30
    epsilon = 1
    delta = 0.1
    r = 0.001

    coher = DPIT.coherence(UT, d)
    L = DPIT.L_nr_it(E[k], E[q1], n)
    noise_variance = DPIT.noise_variance_pow_it(epsilon, p, L, delta)

    eu = DPIT.e_upper(E[k], E[q1], r, d)


    bound = DPIT.epsilon_bound(coher, d, L, E[k], E[q1], epsilon, delta, p)

    eigenvectors, noise = DPIT.power_method(cov, noise_variance, int(L), q1, True)

    res, header = DPIT.generate_result_str(inputfile, epsilon, delta, n, d, q1, p, r, coher, L, noise_variance, bound, noise)
    DPIT.write_summary(res, header, outfile)

    eigenvals = []
    for v in range(eigenvectors.shape[1]) :
        eigenvals.append(DPIT.eigenvalue(cov, eigenvectors[:, v]))

    pd.DataFrame(eigenvals).to_csv(path.dirname(inputfile) + '/eigenvalues_noisy_pit.tsv', sep='\t')
    pd.DataFrame(eigenvectors).to_csv(path.dirname(inputfile) + '/eigenvectors_noisy_pit.tsv', sep='\t')


    #
    eigenvectors, noise = DPIT.power_method(cov, noise_variance, int(L), q1, False)
    eigenvals = []
    for v in range(eigenvectors.shape[1]) :
        eigenvals.append(DPIT.eigenvalue(cov, eigenvectors[:, v]))

    pd.DataFrame(eigenvals).to_csv(path.dirname(inputfile) + '/eigenvalues_pit.tsv', sep='\t')
    pd.DataFrame(eigenvectors).to_csv(path.dirname(inputfile) + '/eigenvectors_pit.tsv', sep='\t')