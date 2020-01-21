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

import convenience as cv
import comparison as co

class AngleRunner():
    def __init__(self, file = '/tmp/'):
        self.ddppca = Distributed_DP_PCA(file = file)
        self.importer = CustomDataImporter()
        self.simulation = SimulationRunner()
        self.outlier = OutlierRemoval()



    def unqeal_split(self, data, interval_end, scale_variance=True, center=True, scale01=False, scale_unit=False, ndims=100, directory='', header=0, rownames=0,  scale_var=True,transpose=False, sep = '\t'):
        """
        This function simulates a multisite PCA with each site having
        varying number of samples.
        :param data:
        :param sites: number of sites to split the data into
        :return:
        """

        print('Shuffling data')
        np.random.shuffle(data)

        Ac = []
        s = 0
        start = 0

        for i in range(len(interval_end)):
            s = s + 1
            end = int(interval_end[i])
            # slice matrix
            data_sub, var_names = self.importer.drop0Columns(data[start:end, :], None, drop=False, noise=True)
            
            print('Local PCA for outlier identification')
            pca, W1, E1 = self.simulation.run_standalone(data_sub, outfile, dims=ndims, header=header, rownames=rownames,center=center,scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.loc', drop_samples=[])

            print('Outlier removal')
            outliers = self.outlier.outlier_removal_mad(pca, 6, 3)
            data_sub = np.delete(data_sub, outliers, 0)


            print('Drop empty variable, scale local data, compute covariance')
            data_sub, vn = self.importer.drop0Columns(data_sub,None, drop = False, noise=True)
            data_sub = self.importer.scale_data(data_sub, center=center, scale_var=scale_variance, scale01=scale01,
                                                scale_unit=scale_unit)
            print('calculating cov')
            noisy_cov = self.ddppca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, nrSamples=data.shape[0],
                                                      nrSites=len(interval_end), noise=False)  # add noise
            start = int(interval_end[i])
            print('finished')


            Ac.append(self.ddppca.perform_SVD(noisy_cov, ndims))



        print('Aggregate local PCAs')
        W, X = self.ddppca.aggregate_partial_SVDs(Ac)
        W = self.ddppca.normalize_eigenvectors(W)
        return (W, X)


    def run_and_compare_unequal(self, datafile, outfile=None, dims=100, header=0, rownames=0, center=True, scale_var=True, scale01=False, scale_unit=False,transpose=False, sep = '\t', reported_angles = 20):

        study_id = path.basename(path.dirname(datafile))
        print('Standalone PCA before outlier removal')
        PCA = SimulationRunner()
        pca, W1, E1 = PCA.run_standalone(datafile, outfile, dims=dims, header=header, rownames=rownames, center=center, scale_var=scale_var, scale01=scale01,scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.before_outlier_removal', log = True)

        print('Logging outliers')
        outliers = self.outlier.outlier_removal_mad(pca, 6, 3)
        with open(outfile + '/removed_outliers.tsv', 'a+') as handle:
            handle.write(cv.collapse_array_to_string(outliers, study_id))

        print('Standalone PCA after outlier removal')
        pca, W1, E1 = PCA.run_standalone(datafile, outfile, dims=dims, header=header, rownames=rownames, center=center,
                                         scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                         transpose=transpose, sep=sep, filename='/pca.after_outlier_removal',
                                         drop_samples=outliers)
        W1 = self.ddppca.normalize_eigenvectors(W1)


        # import data
        data, varnames, sampleids = self.importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile, transpose=False, sep=sep)
        n = data.shape[0]

        interval_end = self.make_test_intervals(n)

        for ar in interval_end:
            for i in range(10): # run it more often the smaller the splits
                print('Current split')
                print(ar)
                W, X = self.unqeal_split(data, ar, ndims=dims)
                angles = self.compute_angles(W1, W, reported_angles)

                #write results
                meta = [len(ar)] + ar
                with open(outfile + '/angles_unequal_splits.tsv', 'a+') as handle:
                    handle.write(cv.collapse_array_to_string(angles, study_id= study_id))
                with open(outfile + '/meta_splits.tsv', 'a+') as handle:
                    handle.write(cv.collapse_array_to_string(meta,  str(i)))
                with open(outfile + '/eigenvalues.tsv', 'a+') as handle:
                    handle.write(cv.collapse_array_to_string(X[1:reported_angles], str(i)))

    def make_test_intervals(self, n):
        # hardcoded for now
        # more extreme cases maybe later
        unequal_splits = list()
        unequal_splits.append([[1.0]])
        unequal_splits.append([[0.1, 0.9], [0.3, 0.7], [0.5, 0.5]])
        unequal_splits.append([[0.2, 0.2, 0.2, 0.2, 0.2], [0.1, 0.1, 0.2, 0.2, 0.4], [0.1, 0.1, 0.1, 0.1, 0.6],
                               [0.2375, 0.2375, 0.2375, 0.2375, 0.05]])

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
        return interval_end

    def compute_angles(self, canonical, split, reported_angles=20):
        angles = list()
        for i in range(min(reported_angles, min(canonical.shape[1], split.shape[1]))):
            angle = co.angle(canonical[:, i], split[:, i])
            angles.append(angle)
        return(angles)


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
    summaryfile = cv.make_eigenvector_path(outfile, path.basename(path.dirname(inputfile)))
    sim.run_and_compare_unequal(inputfile, summaryfile, dims = 20, scale_unit=False, sep = '\t')

