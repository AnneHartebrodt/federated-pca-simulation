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

class Distributed_DP_PCA():

    def __init__(self, logger = None):
        if logger is None:
            self.logger = logging.getLogger('default')
        else:
            self.logger = logger

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
        V,X,W = sc.linalg.svd(Ac)
        W = np.transpose(W)
        W = self.normalize_eigenvectors(W)
        self.logger.info('...done')
        return W[:, 0:ndims],X[0:ndims]

    def local_PCA(self, original, epsilon0, delta0, noise=True, ndims = 10):
        noisy_cov = self.compute_noisy_cov(original, epsilon0, delta0, noise=noise)
        PC = self.perform_SVD(noisy_cov, ndims)
        return PC

    def simulate_multisite_PCA(self, data, sites, epsilon=0.01, delta=0.01, noise=True, ndims=4, scale=True,center=True):
        """
        This function simulates a multisite PCA with each site having
        (almost) the same number of samples.
        :param data:
        :param sites: number of sites to split the data into
        :return:
        """
        print('Running simulation with noise> '+str(noise))
        Ac = []
        s = 0
        start = 0
        interval = int(np.floor(data.shape[0] / sites))
        #data = scale_data(data, center=True, scale_variance=True, scale01=False)
        for i in range(sites):
            s=s+1
            end = min(start + interval, data.shape[0])
            # slice matrix
            data_sub = self.scale_data(data[start:end, :], center=center, scale_variance=scale)
            so = data_sub.shape[0]
            noisy_cov = self.compute_noisy_cov(data_sub, epsilon0= epsilon, delta0=delta, nrSamples=data.shape[0], nrSites=sites, noise=noise) #add noise
            start = start + interval
            Ac.append(self.perform_SVD(noisy_cov, ndims))
        #print(Ac)
        W,X = self.aggregate_partial_SVDs(Ac, ndims=ndims)
        W = self.normalize_eigenvectors(W)

        return(W,X)

    def run_distributed_PCA_locally(self, datasets, epsilon=0.01, delta=0.01, noise=True, ndims=4, scale=True,center=True):
        '''
        This function takes a list of datasets and runs a distributed PCA
        :return:
        '''
        Ac = []
        for data_sub in datasets:
            data_sub = self.scale_data(data_sub, center=center, scale_variance=scale)
            noisy_cov = self.compute_noisy_cov(data_sub, epsilon0= epsilon, delta0=delta, nrSamples=data.shape[0], nrSites=sites, noise=noise) #add noise
            Ac.append(self.perform_SVD(noisy_cov, ndims))
        #print(Ac)
        W,X = self.aggregate_partial_SVDs(Ac, ndims=ndims)
        W = self.normalize_eigenvectors(W)
        return(W,X)

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

    def scale_data(self, data, center=True, scale_variance=False, scale01=False):
        """
        This function centers the data by subtracting the column menans. Scaling to equal variance
        is done by divinding the entries by the column standard deviation. Scaling to values between 0 and 1
        has to be done, in order to assure, that the privacy garantuee is met.
        :param data: nxd numpy array , containing n samples with d variables
        :param center: if true, data is centered by subtracting the column mean
        :param scale_variance: if true, data is scaled by dividing by the standard deviation
        :param scale01: if true, columns are scaled to contain values between 0 and 1
        :return: the scaled data
        """
        #data.dtype = 'float'
        data = data.astype('float')

        if scale_variance or center:
            for column in range(data.shape[1]):
                mean_val = sc.mean(data[:,column])
                var_val = sc.var(data[:,column])
                for elem in range(len(data[:,column])):
                    if center:
                        data[elem, column] = data[elem, column]- mean_val
                    if scale_variance:
                        data[elem, column] = data[elem, column]/var_val
        if scale01:
            max_val = data.max()
            min_val = data.min()
            if(math.isclose(max_val,min_val, rel_tol=1E-15)):
                raise ZeroDivisionError('Insufficient variance in data matrix.')
            for column in range(data.shape[1]):
                for elem in range(len(data[:, column])):
                    data[elem, column] = (data[elem, column] - min_val) / (max_val - min_val)
        return data

    def data_import(self, filename, seed=11, nr_samples=None, header=None, rownames=None, sep='\t', outfile=None):
        """
        Import data, drop 0 variance columns, scale data
        :param filename:
        :param ddpca:
        :param header: Header column? Same as in pandas read csv: None if no header row, row number of there is a header row
        :param rownames: specify the sample id colum: None, if no ids, column number if there are ids
        :return:
        """
        self.logger.info('Reading data (master) ...')
        data = pd.read_csv(filepath_or_buffer=filename, header=header, sep=sep)
        if rownames is not None:
            sample_ids = data.iloc[:, rownames]
            if header is not None:
                data = data.drop(data.columns[rownames], axis=1)
            else:
                data = data.drop([rownames], axis=1)
        else:
            sample_ids = None
        variable_names = data.columns.values
        data = data.values
        data, variable_names = self.drop0Columns(data, variable_names)
        data = self.scale_data(data, True, True, True)
        if nr_samples is not None:
            rnd.seed(seed)
            rnd.choice(data.shape[1], nr_samples)
            print('sample')
        sample_ids = []
        if outfile is not None:
            try:
                self.logger.info('Creating file...')
                os.makedirs(outfile, exist_ok=True)
                self.logger.info('...done')
            except OSError:
                print("Creation of the directory %s failed" % outfile)
            # Transform data to pandas dataframe and print them to the specified location
            pd.DataFrame(data).to_csv(path_or_buf=outfile + '/scaled.data.tsv', sep='\t', header=False, index=False)
            pd.DataFrame(sample_ids).to_csv(path_or_buf=outfile + '/sample.ids.tsv', sep='\t', header=False,
                                            index=False)
            pd.DataFrame(variable_names).to_csv(path_or_buf=outfile + '/variable.names.tsv', sep='\t', header=False,
                                                index=False)
        return (data, sample_ids, variable_names)

    def test_dummy(self):
        self.logger.info('Executing dummy')

    def drop0Columns(self, data, variable_names):
        '''
        Remove columns that have a 0 mean or a zero variance.
        :param data: Data table
        :return: data without 0 columns
        '''

        indices = []
        for column in range(data.shape[1]):
            mean_val = sc.mean(data[:, column])
            var_val = sc.std(data[:, column])
            if (math.isclose(mean_val, 0, rel_tol=1E-15) | math.isclose(var_val, 0, rel_tol=1E-15)):
                indices.append(column)
        data = sc.delete(data, indices, 1)
        variable_names = sc.array(variable_names)
        variable_names = sc.delete(variable_names, indices)
        return (data, variable_names)


if __name__=="__main__":

    filelist = []
    dppca = Distributed_DP_PCA()
    datasets=[]
    for file in filelist:
        data = dppca.data_import(file)
        datasets.append(data)
        dppca.run_distributed_PCA_locally(datasets)

