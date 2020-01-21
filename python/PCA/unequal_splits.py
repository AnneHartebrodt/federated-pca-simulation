import numpy as np
import pandas as pd
import scipy as sc
import random as r
import math as math
from master import Distributed_DP_PCA
from PCA_runner import SimulationRunner
from outlier_removal import OutlierRemoval
import argparse as ap
import numpy as np
import scipy as sc
import copy as copy
import os as os
import pandas as pd
import importlib.util
import scipy.linalg as la
import math as math
from import_export.import_data import CustomDataImporter
import os.path as path
import time as time

import scipy.sparse.linalg as lsa
import numpy as np
import scipy.linalg as la

class AngleRunner():
    def __init__(self, file = '/tmp/'):
        self.ddppca = Distributed_DP_PCA(file = file)
        self.importer = CustomDataImporter()

    def angle(self, v1, v2):
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
        angle = math.degrees(theta)
        # angle can be Nan
        # print(angle)
        if math.isnan(angle):
            return angle
        # return the canonical angle
        if angle > 90:
            return np.abs(angle - 180)
        else:
            return angle


    def unqeal_split(self, data, interval_end, scale_variance=True,
                     center=True, scale01=False, scale_unit=False, ndims=20, directory='', header=0, rownames=0,  scale_var=True,transpose=False, sep = '\t'):
        """
        This function simulates a multisite PCA with each site having
        varying number of samples.
        :param data:
        :param sites: number of sites to split the data into
        :return:
        """
        print('unequal splits')
        PCA = SimulationRunner()
        OR = OutlierRemoval()

        np.random.shuffle(data)

        Ac = []
        s = 0
        start = 0
        # data = scale_data(data, center=True, scale_variance=True, scale01=False)
        for i in range(len(interval_end)):
            s = s + 1
            end = int(interval_end[i])
            # slice matrix
            data_sub, var_names = self.importer.drop0Columns(data[start:end, :], None, drop=False, noise=True)
            
            print('regular pca')
            pca, W1, E1 = PCA.run_standalone(data_sub, outfile, dims=ndims, header=header, rownames=rownames,
                                             center=center,
                                             scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                             transpose=transpose, sep=sep, filename='/pca.loc',
                                             drop_samples=[])
            
	    
            outliers = OR.outlier_removal_mad(pca, 6, 3)
            print('outliers')
            print(outliers)
            data_sub = np.delete(data_sub, outliers, 0)

            data_sub, vn = self.importer.drop0Columns(data_sub,None, drop = False, noise=True)
            data_sub = self.importer.scale_data(data_sub, center=center, scale_var=scale_variance, scale01=scale01,
                                                scale_unit=scale_unit)
            print('calculating cov')
            noisy_cov = self.ddppca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, nrSamples=data.shape[0],
                                                      nrSites=len(interval_end), noise=False)  # add noise
            start = int(interval_end[i])
            print('finished')


            Ac.append(self.ddppca.perform_SVD(noisy_cov, ndims))
            # print(Ac)
        print('aggregating')
        W, X = self.ddppca.aggregate_partial_SVDs(Ac, ndims=ndims)
        W = self.ddppca.normalize_eigenvectors(W)

        print('....finished')


        #
        #self.save_PCA(None, W, X, directory + '/pca')
        return (W, X)


    def run_and_compare_unequal(self, datafile, outfile=None, dims=5, header=0, rownames=0, center=True, scale_var=True, scale01=False, scale_unit=False,transpose=False, sep = '\t'):

        study_id = path.basename(path.dirname(datafile))
        print('standalone before outlier removal')
        PCA = SimulationRunner()
        pca, W1, E1 = PCA.run_standalone(datafile, outfile, dims=dims, header=header, rownames=rownames, center=center, scale_var=scale_var, scale01=scale01,scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.before_outlier_removal')

        ore = OutlierRemoval()
        outliers = ore.outlier_removal_mad(pca, 6, 3)
        print(outliers)
        print('standalone after outlier removal')
        print(datafile)
        pca, W1, E1 = PCA.run_standalone(datafile, outfile, dims=dims, header=header, rownames=rownames, center=center,
                                         scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                         transpose=transpose, sep=sep, filename='/pca.after_outlier_removal',
                                         drop_samples=outliers)


        W1 = self.ddppca.normalize_eigenvectors(W1)

        data, varnames, sampleids = self.importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile, transpose=False, sep=sep)

        n = data.shape[0]

        # hardcoded for now
        # more extreme cases maybe later
        unequal_splits = list()
        unequal_splits.append([[1.0]])
        unequal_splits.append([[0.2, 0.8], [0.3, 0.7], [0.4, 0.6]])
        unequal_splits.append([[0.2, 0.4, 0.4], [0.2, 0.2, 0.6], [0.15, 0.35, 0.5]])
        unequal_splits.append([[0.2, 0.2, 0.2, 0.4], [0.15, 0.15, 0.2, 0.5], [0.1, 0.25, 0.25, 0.4]])

        interval_end = list()
        sum = 0
        for i in unequal_splits:
            for j in i:
                inter = list()
                for k in j:
                    sum = sum + k
                    inter.append(np.ceil(n * sum))
                sum = 0
                interval_end.append(inter)

        for ar in interval_end:
            for i in range(len(ar)+1): # run it more often the smaller the splits
                print(ar)
                W, X = self.unqeal_split(data, ar)
                angles = self.compute_angles(W1, W)
                angles_str = self.collapse_array_to_string(angles, study_id= study_id)
                meta = [len(ar)] + ar

                meta = self.collapse_array_to_string(meta,  str(i))
                eig_str = self.collapse_array_to_string(X[1:20], str(i))
                with open(outfile + '/angles_unequal_splits.tsv', 'a+') as handle:
                    handle.write(angles_str)
                with open(outfile + '/meta_splits.tsv', 'a+') as handle:
                    handle.write(meta)
                with open(outfile + '/eigenvalues.tsv', 'a+') as handle:
                    handle.write(eig_str)





    def compute_angles(self, canonical, split):
        angles = list()
        for i in range(min(20, min(canonical.shape[1], split.shape[1]))):
            angle = self.angle(canonical[:, i], split[:, i])
            angles.append(angle)
        return(angles)




    def extract_eigenvals(self,E):
        '''
        Eigendecomposition from scipy.linalg.sparse returns eigenvalues ordered in
        increasing order, followed by eigenvalues which are 0.
        Eigenvalues are returned in decreasing order ommiting th 0s alltogether
        Args:
            E: Eigenvalues from a sparse singular value decomposition.

        Returns: Eigenvalue vector in decreasing order, without 0s.

        '''
        E = E[E != 0]
        return E, indz

    def collapse_array_to_string(self, a, study_id):
        res = study_id  + '\t'
        for e in a:
            res = res + str(e) + '\t'
        res = res + '\n'
        return res

    def make_eigenvector_path(self, inputfile, foldername):
        """
        creates a folder called eigenvectors in the input directory
        if the eigenvector folder already exists a folder named
        'eigenvectors<currentTimestamp> will be created
        Args:
            inputfile: name of the inputfile

        Returns: a pathname according to the stated rules

        """
        print(path.dirname(inputfile))
        if not os.path.exists(path.dirname(inputfile)+'/'+foldername):
            pn = path.dirname(inputfile)+'/'+foldername
            os.makedirs(pn)
        else:
            print('Path exists')
            pn = path.dirname(inputfile)+'/'+foldername+str(time.process_time())
            os.makedirs(pn)
        return pn


if __name__=="__main__":
    print('run split script')

    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    args = parser.parse_args()

    inputfile = args.f
    outfile = args.o

    #inputfile ='/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_trunc.tsv'
    #outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/testttt/'

    cd = CustomDataImporter()
    sim = AngleRunner()
    summaryfile = sim.make_eigenvector_path(outfile, path.basename(path.dirname(inputfile)))
    sim.run_and_compare_unequal(inputfile, summaryfile, dims = 20, scale_unit=False, sep = '\t')

