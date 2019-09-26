from PCA.master import Distributed_DP_PCA
import argparse as ap
import numpy as np
import scipy as sc
import copy as copy
import os as os
import pandas as pd
import importlib.util



class SimulationRunner():
    def __init__(self):
        spec = importlib.util.spec_from_file_location("module.name","../../import_export/import_data.py")
        foo = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(foo)
        self.ddppca = Distributed_DP_PCA()
        self.importer = foo.CustomDataImporter()


    def create_annotation(self, n, split, i, epsilon, delta):
        proj = np.ones((n, 1)) * split
        proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * i)), axis=1)
        proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * epsilon)), axis=1)
        proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * delta)), axis=1)
        return proj


    def run_multiple_simulations(self, datafile, dims, noise, splits, repeat, epsilons=[], deltas=[],
                                 dirname='', save_eigen=False, center = True, scale_var = True, scale01=False,
                                 scale_unit=True, transpose=False, header=None, rownames=None):

        # import data
        data, sample_ids, variable_names = self.importer.data_import(datafile, header, rownames,outfile=dirname, transpose=transpose)
        scaled_data= self.importer.scale_data(copy.deepcopy(data), center=center, scale_unit=scale_unit, scale_var=scale_var, scale01=scale01)
        # make result dir
        try:
            os.mkdir(dirname)
        except OSError:
            print("Creation of the directory %s failed" % dirname)


        # fill a dummy value for noiseless simulation
        if not noise:
            epsilons = [2]
            deltas = [2]
            repeat = 1
        # for a noisy simulation generate some example epsilons
        if noise and len(epsilons)==0:
            epsilons = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 100]
        if noise and len(deltas)==0:
            deltas = [0.01]

        print('Running simulation with noise: ' + str(noise))
        backup = copy.deepcopy(data)
        dm = min(dims, data.shape[1])
        results = np.empty(shape=(1, dm + 4))
        eigenvalues = np.empty(shape=(1, dm + 4))
        eigenvectors = np.empty(shape=(1, dm + 4))
        for epsilon in epsilons:
            for delta in deltas:
                v1 = None
                for i in range(repeat):
                    for split in range(1, splits + 1):

                        data = copy.deepcopy(backup)
                        vec, val = self.simulate_multisite_PCA(data, split, epsilon, delta, noise=noise, ndims=dims,
                                                                scale_variance=scale_var, center=center, scale01=scale01,
                                                               scale_unit=scale_unit)
                        if v1 is None:
                            v1 = copy.deepcopy(vec)
                        else:
                            vec = self.ddppca.normalize_eigenspaces([v1, vec])[1]

                        # projection matrix for sanity check.
                        proj = sc.dot(scaled_data, vec[:, 0:dm])
                        proj = np.concatenate((proj, self.create_annotation(proj.shape[0], split, i, epsilon, delta)),
                                              axis=1)
                        results = np.concatenate((results, proj), axis=0)
                        # eigenvalues
                        ar = np.array(np.concatenate((val, np.array([split, i, epsilon, delta])))).reshape((1, dm + 4))
                        eigenvalues = np.concatenate((eigenvalues, ar), axis=0)
                        # eigenvectors
                        vec = np.concatenate((vec, self.create_annotation(vec.shape[0], split, i, epsilon, delta)), axis=1)
                        eigenvectors = np.concatenate((eigenvectors, vec), axis=0)

        # remove first, random row
        results = np.delete(results, 0, 0)
        eigenvectors = np.delete(eigenvectors, 0, 0)
        eigenvalues = np.delete(eigenvalues, 0, 0)

        if noise:
            filename = dirname +'noise_pca'
        else:
            filename = dirname +'no_noise_pca'

        self.save_PCA(results, eigenvectors, eigenvalues, filename)
        return results, eigenvalues, eigenvectors


    def simulate_multisite_PCA(self, data, sites, epsilon=0.01, delta=0.01, noise=True, ndims=4, scale_variance=True,
                               center=True, scale01 = False, scale_unit=True):
        """
        This function simulates a multisite PCA with each site having
        (almost) the same number of samples.
        :param data:
        :param sites: number of sites to split the data into
        :return:
        """
        print('Running simulation with noise> ' + str(noise))
        Ac = []
        s = 0
        start = 0
        interval = int(np.floor(data.shape[0] / sites))
        # data = scale_data(data, center=True, scale_variance=True, scale01=False)
        for i in range(sites):
            s = s + 1
            end = min(start + interval, data.shape[0])
            # slice matrix
            data_sub = self.importer.scale_data(data[start:end, :], center=center, scale_var=scale_variance, scale01=scale01, scale_unit=scale_unit)
            so = data_sub.shape[0]
            noisy_cov = self.ddppca.compute_noisy_cov(data_sub, epsilon0=epsilon, delta0=delta, nrSamples=data.shape[0],
                                               nrSites=sites, noise=noise)  # add noise
            start = start + interval
            Ac.append(self.ddppca.perform_SVD(noisy_cov, ndims))
        # print(Ac)
        W, X = self.ddppca.aggregate_partial_SVDs(Ac, ndims=ndims)
        W = self.ddppca.normalize_eigenvectors(W)

        return (W, X)


    def run_distributed_PCA_locally(self, datasets, epsilon=0.01, delta=0.01, noise=True, ndims=4, scale=True,
                                    center=True, sites=1):
        '''
        This function takes a list of datasets and runs a distributed PCA
        :return:
        '''
        Ac = []
        for data_sub in datasets:
            data_sub = self.importer.scale_data(data_sub, center=center, scale_variance=scale)
            noisy_cov = self.ddppca.compute_noisy_cov(data_sub, epsilon0=epsilon, delta0=delta, nrSamples=data.shape[0],
                                               nrSites=sites, noise=noise)  # add noise
            Ac.append(self.ddppca.perform_SVD(noisy_cov, ndims))
        # print(Ac)
        W, X = self.ddppca.aggregate_partial_SVDs(Ac, ndims=ndims)
        W = self.ddppca.normalize_eigenvectors(W)
        return (W, X)

    def run_standalone(self, datafile, outfile=None, dims=5, header=None, rownames=4, center=True, scale_var=True, scale01=False, scale_unit=True,
                       transpose = False):
            '''
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
            '''
            data, sample_ids, variable_names = self.importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile, sep='\t', transpose=transpose)
            data = self.importer.scale_data(data, center, scale_var, scale01, scale_unit)
            # standalone PCA
            pca, W, s = self.ddppca.standalone_pca(data, ndims=dims)
            self.save_PCA(pca,W,s, outfile+'/pca')
            return pca, W, s

    def save_PCA(self, pca, W, s, outfile):
        pd.DataFrame(pca).to_csv(outfile+'.projection', sep='\t', header=None, index=False)
        pd.DataFrame(W).to_csv(outfile + '.eigenvectors', sep='\t', header=None, index=False)
        pd.DataFrame(s).to_csv(outfile + '.eigenvalues', sep='\t', header=None, index=False)


if __name__=="__main__":
    print('run')

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
    parser.add_argument('-v',action='store_true',
                        help='Scale variables to unit variance', default=True)
    parser.add_argument('-u',action='store_true',
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
        deltas = parse_array(args.g)
    except:
        print('Incompatible epsilon or delta parameters')
        print('default to epsilon=0.1, delta = 0.001')
        epsilons = [0.1]
        deltas = [0.001]

    print('file: ' + str(args.f))
    print('dimensions: ' + str(args.d))
    print('output directory: ' + str(args.p))
    print('number of hospitals: ' + str(args.k))
    print('save eigenvalues: ' + str(args.s))
    print('repeats: ' + str(args.r))
    print('column headers: ' + str(args.c))
    print('sample ids: ' + str(args.i))

    #simulation.run_standalone(datafile, outfile=outfile, dims=dimensions, header=header,rownames=rownames,center = center, scale_var = scale_var, scale01 = scale01, scale_unit=scale_unit, transpose = True)
    #simulation.run_multiple_simulations(datafile=datafile, dims=dimensions, noise=noise, splits=splits, repeat = repeats, epsilons=[0.01], deltas=[0.01],dirname=outfile, save_eigen=False, transpose = True, center = center, scale_var = scale_var, scale01 = scale01, scale_unit=scale_unit)

    simulation = SimulationRunner()

    if args.A:
        simulation.run_standalone(args.f, outfile=args.p, dims=args.d, header=header, rownames=rownames,
                              center=args.t, scale_var=args.v, scale01=args.z, scale_unit=args.u,
                              transpose=False)
    if args.B:
        simulation.run_multiple_simulations(datafile=args.f, dims=args.d, repeat=args.r, header=header, rownames=rownames,
                                        epsilons=epsilons, deltas=deltas, dirname=args.p, save_eigen=args.s,
                                        transpose=False, center=args.t, scale_var=args.v, scale01=args.z,
                                        scale_unit=args.u, noise=False, splits=args.k)

    if args.C:
        simulation.run_multiple_simulations(datafile=args.f, dims=args.d, repeat=args.r, header=header, rownames=rownames,
                                        epsilons=epsilons, deltas=deltas, dirname=args.p, save_eigen=args.s,
                                        transpose=False, center=args.t, scale_var=args.v, scale01=args.z,
                                        scale_unit=args.u, noise=True, splits=args.k)



