import os.path as path

import scipy.linalg as lsa

#import python.PCA.import_export.
import python.PCA.horizontal.power_iteration as powerit
import python.PCA.horizontal.proxy_covariance as dpca
import python.PCA.vertical.vertical_pca_library as vert
# import import_export.easy_import as easy
import argparse as ap
import numpy as np
import scipy.sparse.linalg as lsa

def run_standalone(data, outfile=None, p=10, header=None, rownames=None, center=True, scale_var=True, scale01=False,
                   scale_unit=True, transpose=False, sep='\t', drop_samples=[], log=True, exp_var=0.5, epsilon=1,
                   delta=1, noise=False):
    """
    This function performs a regular principal component analysis and saves the result to files containing
    the projection the
    :param datafile: Unscaled datafile
    :param ddpca:
    :param outfile: path and name of the output file without extension
    :param dims: Number of dimensions to return (#eigenvectors and corresponding eigenvalues)
    :param seed: random seed
    :param nr_samples: #variables to select if not all the data columns are to be used for the pca
    :param header: row number which contains the header/ number of header rows
    :param rownames: column number which contains the rownames/sample ids
    :return: projection, eigenvectors and eigenvalues
    """
    # if data is a string, data needs to be read first, otherwise it is
    # assumed to be scaled and ready to use
    if isinstance(data, str):
        data = easy.easy_import(data, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,
                                scale01=scale01, scale_unit=scale_unit, drop_samples=drop_samples, log=log,
                                outfile=outfile, transpose=transpose)
    # run pca
    noisy_cov = dpca.compute_noisy_cov(data, epsilon0=1, delta0=1, noise=False)
    if noise:
        params = powerit.determine_parameters(data, epsilon, delta)
        pca, W, s, nriter, noise = powerit.power_method(noisy_cov, params['sigma'], params['L'], params['p'], noise)
    else:
        noise_variance = 1
        params = []
        pca, W, s, nr_iter, noise = powerit.power_method(noisy_cov, sigma=noise_variance, p=p, noise=noise)
    # return the projected datapoints, the eigenvectors and the eigenvalues
    return pca, W, s, nr_iter, params


def simulate_distributed(local, p=8, weights=None, tolerance=0.000001, scipy = None):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[1]
    X = powerit.generate_random_gaussian(d, p, sigma=1)
    X, R = lsa.qr(X, mode='economic')
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    while not converged and count<1000:
        count = count + 1
        locals = []
        for l in local:
            d1 = np.dot(l, X)
            d1 = np.dot(l.T, d1)
            locals.append(d1)
        X, E, converged = powerit.pooling_step(locals, X, weights=weights)
    ord = np.argsort(E)
    X = np.flip(X[:, ord], axis=1)
    E = np.flip(np.sort(E))
    return X, E, count


def simulate_distributed_multi_local(local, p=8, weights=None, nr_local_rounds = 5, tolerance=0.000001):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[0]
    X = powerit.generate_random_gaussian(d, p, sigma=1)
    X, R = lsa.qr(X, mode='economic')
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    X_prev = X
    while not converged:
        count = count + 1
        locals = []
        for l in local:
            # run a few iterations locally before sending the result to the aggrgator
            for i in range(nr_local_rounds):
                locals.append(powerit.local_step(X, l))
        X, E, converged = powerit.pooling_step(locals, X, weights=weights)

        converged, conv, converged_eigenvals, delta = vert.convergence_checker(X, X_prev)
        print(converged)

    ord = np.argsort(E)
    X = np.flip(X[:, ord], axis=1)
    E = np.flip(np.sort(E))
    return X, E, count


def simulate_deflation(data, outfile=None, dims=10, header=None, rownames=None, center=True, scale_var=True,
                       scale01=False, scale_unit=True, transpose=False, sep='\t', filename='/pca', drop_samples=[],
                       log=True, exp_var=0.5, epsilon=1, delta=1, noise=False):
    """
    This function performs a regular principal component analysis and saves the result to files containing
    the projection the
    :param datafile: Unscaled datafile
    :param ddpca:
    :param outfile: path and name of the output file without extension
    :param dims: Number of dimensions to return (#eigenvectors and corresponding eigenvalues)
    :param seed: random seed
    :param nr_samples: #variables to select if not all the data columns are to be used for the pca
    :param header: row number which contains the header/ number of header rows
    :param rownames: column number which contains the rownames/sample ids
    :return: projection, eigenvectors and eigenvalues
    """

    # if data is a string, data needs to be read first, otherwise it is
    # assumed to be scaled and ready to use
    if isinstance(data, str):
        data = easy.easy_import(data, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,
                                scale01=scale01, scale_unit=scale_unit, drop_samples=drop_samples, log=log,
                                outfile=outfile, transpose=transpose)
    # run pca
    noisy_cov = dpca.compute_noisy_cov(data, epsilon0=epsilon, delta0=delta, noise=False)
    eigenvalues = []
    eigenvectors = []
    nr_iter = []
    for d in range(dims):

        pca, W, E, nriter= powerit.power_iteration(noisy_cov)
        noisy_cov = powerit.hotelling_deflation(noisy_cov, W, E, False)
        eigenvalues.append(E)
        eigenvectors.append(W.flatten())
        nr_iter.append(nriter)
        print(d)

    eigenvectors = np.column_stack(eigenvectors)

    return eigenvectors, eigenvalues, nr_iter


def generate_result_str(inputfile, epsilon, delta, n, d, q, p, r, coherence, L, noise_variance, e, noise_vals,
                        assumed_noise):
    res = path.basename(path.dirname(inputfile)) + '\t' + str(epsilon) + '\t' + str(delta) + '\t' + str(n) + '\t' + str(
        d) + '\t' + str(q) + '\t' + str(p) \
          + '\t' + str(r) + '\t' + str(coherence) + '\t' + str(L) + '\t' + str(noise_variance) + '\t' + str(
        e) + '\t' + str(assumed_noise) + '\t'
    header = 'study.id\tepsilon\tdelta\tnr.samples\tnr.dimension\tnr.itermediate\tit.rank\tr\tcoherence\tnr.iterations\tnoise.variance\terror.bound\tassumed.noise\t'
    i = 1
    for v in noise_vals:
        res = res + str(v) + '\t'
        header = header + 'n1.' + str(i) + '\t'
        i = i + 1
    return (res, header)


if __name__ == '__main__':
    print('power iteration runner')
    parser = ap.ArgumentParser(description='Run distributed PCA simulation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-d', metavar='dimensions', type=int, help='number of principal components to return')
    parser.add_argument('-p', metavar='output directory', type=str, help='output directory for simulation study.')
    parser.add_argument('-k', metavar='number_hospitals', type=int, help='Number of simulated hospitals', default=5)
    parser.add_argument('-s', action='store_true',
                        help='If true the generated eigenspaces are saved (!a lot of large files!)', default=False)
    parser.add_argument('-r', metavar='Repeats', type=int, help='Number of repetitions of the sampling process',
                        default=5)
    parser.add_argument('-c', help='True of data has column headers',
                        default=False, action='store_true')
    parser.add_argument('-i', metavar='sampleids', type=int,
                        help='Dataframe column which contains the sample ids, 0 index,', default=-1)
    parser.add_argument('-e', metavar='epsilons', type=str, help='epsilons to simulate separated by a comma')
    parser.add_argument('-g', metavar='deltas', type=str, help='deltas to simulate separated by a comma')
    parser.add_argument('-v', action='store_true',
                        help='Scale variables to unit variance', default=True)
    parser.add_argument('-u', action='store_true',
                        help='Scale samples to unit norm', default=True)
    parser.add_argument('-z', action='store_true',
                        help='Scale variables between 0 and 1', default=False)
    parser.add_argument('-t', action='store_true',
                        help='Center variables by substracting the mean', default=True)

    parser.add_argument('-A', action='store_true',
                        help='Run standalone simulation', default=False)
    parser.add_argument('-B', action='store_true',
                        help='Run distributed simulation without noise', default=False)
    parser.add_argument('-C', action='store_true',
                        help='Run distributed simulation with noise', default=False)
    parser.add_argument('-D', action='store_true',
                        help='Run distributed, submit a list of files (comma sep), noise = F', default=False)
    parser.add_argument('-E', action='store_true',
                        help='Run distributed, submit a list of files (comma sep), noise = T', default=False)
    parser.add_argument('-V', action='store_true',
                        help='Save Covariance matrix before and after noise', default=False)
    args = parser.parse_args()


    def parse_array(value_str):
        values_str = value_str.split(',')
        values_int = []
        for v in values_str:
            values_int.append(float(v))
        return values_int


    if args.c:
        header = 0
    else:
        header = None

    if args.i == -1:
        rownames = None
    else:
        rownames = args.i

    try:
        epsilons = parse_array(args.e)
    except:
        print('Incompatible epsilon parameters')
        print('default to epsilon=0.1')
        epsilons = [0.1]

    try:
        deltas = parse_array(args.g)
    except:
        print('Incompatible delta parameters')
        print('default to delta=0.001')
        deltas = [0.001]

    print('file: ' + str(args.f))
    print('dimensions: ' + str(args.d))
    print('output directory: ' + str(args.p))
    print('number of hospitals: ' + str(args.k))
    print('save eigenvalues: ' + str(args.s))
    print('repeats: ' + str(args.r))
    print('column headers: ' + str(args.c))
    print('sample ids: ' + str(args.i))

    if args.A:
        run_standalone(args.f, outfile=args.p, dims=args.d, header=header, rownames=rownames,
                       center=args.t, scale_var=args.v, scale01=args.z, scale_unit=args.u,
                       transpose=False)
