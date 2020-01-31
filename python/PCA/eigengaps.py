import sys
from import_export.import_data import CustomDataImporter
from PCA_runner import SimulationRunner
import scipy.sparse.linalg as lsa
import numpy as np
import scipy.linalg as la
import os.path as path
import os as os
import scipy as sc
from outlier_removal import OutlierRemoval
import argparse as ap
import pandas as pd



def sigma_analyze_gauss(epsilon, n, delta=None):
    '''

    Args:
        epsilon: Epsilon parameter, multiply by n if covariance matrix is scaled by n
        n: nr of samples (rows) in dataframe
        delta: (optional) if not given, delta will be set as 1/(n*100) to make it cryptogrphically small

    Returns: the variance of the gaussian noise to be added to covariance matrix

    '''
    if delta is None:
        # Assuming d=1/q(N*100)
        sigma = (1 / epsilon) * np.sqrt(2 * np.log(1.25 * 100 * n))
    else:
        sigma = (1 /epsilon) * np.sqrt(2 * np.log(1.25/delta))
    return sigma


def required_eigengap(sigma, n):
    '''

    Args:
        sigma: variance of gaussian noise from Analyze Gauss
        n: number of samples

    Returns: required eigengap

    '''
    eigengap = np.sqrt(n)*sigma
    return eigengap


def generate_result_eigengap(epsilon, E, nr_samples):
    res = path.basename(path.dirname(inputfile)) + '\t' + str(epsilon) + '\t' + str(nr_samples) + '\t'
    header = 'study.id\tepsilon\tnr.samples\t'
    max_index=-1
    E_gap=list()
    for i in range(0, (len(E)-1)):
        E_gap.append(E[i]-E[i+1])
        res+=str(E[i]-E[i+1])+'\t'
        header+='eigengap.'+str(i+1)+'\t'
        if(E[i]-E[i+1]<eigengap and max_index==-2):
            max_index = i-1
    return (res, header)

def write_summary(res, header, outfile):
    try:
        os.makedirs(path.dirname(outfile))
    except OSError:
        print(path.dirname(outfile)+' could not be created')
    else:
        print(path.dirname(outfile)+ ' was created')

    if not (path.exists(outfile)):
        with open(outfile, 'w') as handle:
            handle.write(header+'\n')

    with open(outfile, 'a+') as handle:
        handle.write(res+'\n')

def extract_eigenvals(E):
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


if __name__=="__main__":
    print('run eigengap script')

    parser = ap.ArgumentParser(description='Eigengaps')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    #parser.add_argument('-v', metavar='explained_var', type=float, help='explained variance')
    parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    #parser.add_argument('-m', metavar='mult_dims_ret', type=str, help='comma separated list of intermediate dimensions', default = 1)
    #parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default = 100)
    args = parser.parse_args()

    inputfile = args.f
    outfile = args.o
    #exp_var = args.v
    #mult_dims_ret = args.m
    sep = args.s
    #dims = args.d

    #inputfile ='/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_trunc.tsv'
    #outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/testttt/sum.txtxt'
    exp_var = 0.8
    #sep = ','
    #dims = 100


    dpca = SimulationRunner()
    pca, W, E = dpca.run_standalone(inputfile, outfile=None, dims=1000, header=0, rownames=None, center=True, scale_var=True, scale01=False, scale_unit=True,
                       transpose = False, sep=sep, filename = '/pca', drop_samples =[], log = True, exp_var=0.99)

    n = pca.shape[0]
    epsilons = [0.01, 0.1, 1, 10]
    deltas = [None, 0.01, 0.1, 1]

    name = path.basename(path.dirname(inputfile))
    result = []
    req_eigengaps= []
    for epsilon in epsilons:
        for delta in deltas:
            eigengap = required_eigengap(sigma_analyze_gauss(n * epsilon, n, delta), n)
            req_eigengaps.append(eigengap)
            result.append([name, epsilon, delta, eigengap])

    df = pd.DataFrame(result)
    df.columns = ['study.id', 'epsilon', 'delta', 'req.eigengap']
    df.to_csv(path.dirname(outfile)+'/required.eigengaps.tsv', sep = '\t', mode='a+', index=False)


    res, header = generate_result_eigengap(epsilon, E[0:20],n)
    write_summary(res, header, outfile)

    orm = OutlierRemoval()
    outliers = orm.outlier_removal_mad(pca, 6, 3)
    pca, W, E = dpca.run_standalone(inputfile, outfile=None, dims=1000, header=0, rownames=None, center=True,
                                    scale_var=True, scale01=False, scale_unit=True,
                                    transpose=False, sep=sep, filename='/pca', drop_samples=outliers, log=True, exp_var=0.99)


    res, header = generate_result_eigengap(epsilon, E[0:20],n)
    write_summary(res, header, path.dirname(outfile)+'/eigengaps_aor.txt')