import argparse as ap
import os.path as path

import numpy as np
import pandas as pd

import proxy_covariance_runner as runner
import convenience as cv
import outlier_removal as orm


def sigma_analyze_gauss(epsilon, n, delta=None):
    """
    Args:
        epsilon: Epsilon parameter, multiply by n if covariance matrix is scaled by n
        n: nr of samples (rows) in dataframe
        delta: (optional) if not given, delta will be set as 1/(n*100) to make it cryptogrphically small

    Returns: the variance of the gaussian noise to be added to covariance matrix

    """

    # sigma = (1 / (n * epsilon0)) * sc.sqrt(2 * np.log(1.25 / delta0))
    if delta is None:
        # Assuming d=1/q(N*100)
        sigma = (1 / epsilon) * np.sqrt(2 * np.log(1.25 * 100 * n))
    else:
        sigma = (1 / epsilon) * np.sqrt(2 * np.log(1.25 / delta))
    return sigma


def required_eigengap(sigma, n):
    """
    Args:
        sigma: variance of gaussian noise from Analyze Gauss
        n: number of samples

    Returns: required eigengap
    """

    eigengap = np.sqrt(n) * sigma
    return eigengap


def generate_result_eigengap(epsilon, eigenvalues, nr_samples):
    res = path.basename(path.dirname(inputfile)) + '\t' + str(epsilon) + '\t' + str(nr_samples) + '\t'
    header = 'study.id\tepsilon\tnr.samples\t'
    eigenvalues_gap = list()
    for i in range(0, (len(eigenvalues) - 1)):
        eigenvalues_gap.append(eigenvalues[i] - eigenvalues[i + 1])
        res += str(eigenvalues[i] - eigenvalues[i + 1]) + '\t'
        header += 'eigengap.' + str(i + 1) + '\t'
    return res, header


if __name__ == "__main__":
    print('run eigengap script')
    # parser = ap.ArgumentParser(description='Eigengaps')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    #
    # parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    # args = parser.parse_args()
    #
    # inputfile = args.f
    # outfile = args.o
    # sep = args.s

    inputfile ='/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_trunc.tsv'
    outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/testttt/sum.txtxt'
    exp_var = 0.8
    sep = ','

    pca, W, E = runner.run_standalone(inputfile, outfile=None, dims=1000, header=0, rownames=None, center=True,
                                      scale_var=True, scale01=False, scale_unit=True,
                                      transpose=False, sep=sep, filename='/pca', drop_samples=[], log=True,
                                      exp_var=0.99)

    n = pca.shape[0]
    epsilon=''
    res, head = generate_result_eigengap(epsilon, E[0:20], n)
    cv.write_summary(res, head, outfile)

    epsilons = [0.01, 0.1, 1, 10]
    deltas = [None, 0.01, 0.1, 1]

    name = path.basename(path.dirname(inputfile))
    result = []
    req_eigengaps = []
    for epsilon in epsilons:
        for delta in deltas:
            eigengap = required_eigengap(sigma_analyze_gauss(n * epsilon, n, delta), n)
            req_eigengaps.append(eigengap)
            result.append([name, epsilon, delta, eigengap])

    df = pd.DataFrame(result)
    df.columns = ['study.id', 'epsilon', 'delta', 'req.eigengap']
    df.to_csv(path.join(path.dirname(outfile), 'required.eigengaps.tsv'), sep='\t', mode='a+', index=False)



    outliers = orm.outlier_removal_mad(pca, 6, 3)
    pca, W, E = runner.run_standalone(inputfile, outfile=None, dims=1000, header=0, rownames=None, center=True,
                                      scale_var=True, scale01=False, scale_unit=True,
                                      transpose=False, sep=sep, filename='/pca', drop_samples=outliers, log=True,
                                      exp_var=0.99)
    epsilon  = ''
    res, head = generate_result_eigengap(epsilon, E[0:20], n)
    cv.write_summary(res, head, path.join(path.dirname(outfile), 'eigengaps_aor.txt'))
