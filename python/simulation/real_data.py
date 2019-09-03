from PCA.master import Distributed_DP_PCA
import scipy as sc
import pandas as pd
import argparse as ap
import PCA.simulation_runner as simul
import math
import os
import numpy.random as rnd

def data_import(filename, ddpca, seed=11, nr_samples=None, header=None, rownames=None, sep = '\t', outfile = None):
    """
    Import data, drop 0 variance columns, scale data
    :param filename:
    :param ddpca:
    :param header: Header column? Same as in pandas read csv: None if no header row, row number of there is a header row
    :param rownames: specify the sample id colum: None, if no ids, column number if there are ids
    :return:
    """
    print(nr_samples)
    print(seed)
    data = pd.read_csv(filepath_or_buffer=filename, header=header, sep =sep)
    if rownames is not None:
        sample_ids = data.iloc[:,rownames]
        if header is not None:
            data = data.drop(data.columns[rownames], axis=1)
        else:
            data = data.drop([rownames], axis=1)
    else:
        sample_ids = None
    variable_names = data.columns.values
    data = data.values
    data, variable_names = drop0Columns(data, variable_names)
    data = ddpca.scale_data(data, True, True, True)
    if nr_samples is not None:
        rnd.seed(seed)
        rnd.choice(data.shape[1], nr_samples)
        print('sample')
    sample_ids= []
    if outfile is not None:
        try:
            os.mkdir(outfile)
        except OSError:
            print("Creation of the directory %s failed" % outfile)
        # Transform data to pandas dataframe and print them to the specified location
        pd.DataFrame(data).to_csv(path_or_buf=outfile+'/scaled.data.tsv', sep = '\t', header=False, index=False)
        pd.DataFrame(sample_ids).to_csv(path_or_buf=outfile+'/sample.ids.tsv', sep = '\t', header=False, index=False)
        pd.DataFrame(variable_names).to_csv(path_or_buf=outfile + '/variable.names.tsv', sep='\t', header=False, index=False)
    return (data, sample_ids, variable_names)

def drop0Columns(data, variable_names):
    '''
    Remove columns that have a 0 mean or a zero variance.
    :param data: Data table
    :return: data without 0 columns
    '''

    indices = []
    for column in range(data.shape[1]):
        mean_val = sc.mean(data[:, column])
        var_val = sc.std(data[:, column])
        if(math.isclose(mean_val, 0, rel_tol=1E-15) | math.isclose(var_val, 0, rel_tol=1E-15) ):
            indices.append(column)
    data = sc.delete(data, indices, 1)
    variable_names = sc.array(variable_names)
    variable_names = sc.delete(variable_names, indices)
    return (data,variable_names)

def run_standalone(datafile,ddpca, outfile=None, dims=5, seed = 11, nr_samples=None, header=None, rownames=4):
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
    data, sample_ids, variable_names = data_import(datafile, ddpca, seed, nr_samples, header, rownames)
    # standalone PCA
    # Data already scaled and cenetered
    pca, W, s = simul.perform_standalone_pca(data, ndims=dims, dppca=ddpca, center=False, scale=False, scale01=False)
    pca1 = pd.DataFrame(pca)
    pca1.to_csv(outfile, sep='\t', header=None, index=False)
    W1 = pd.DataFrame(W)
    W1.to_csv(outfile+'.eigenvectors', sep='\t', header=None, index=False)
    s1 = pd.DataFrame(s)
    s1.to_csv(outfile+'.eigenvalues', sep='\t', header = None, index = False)
    return pca, W, s


def run_1_n_sites_with_noise(datafile, ddpca, dirname, noise=True, ndims=4, maxSplits=5, repeat= 10, epsilons=None, deltas=None, save_eigen=False, nr_variables=None, seed = 11, header=None, rownames=None):

    print('file: '+str(datafile))
    print('output file: ')
    print('dimensions: '+str(ndims))
    print('output directory: '+str(dirname))
    print('random seed: '+str(seed))
    print('number of sampled variables: '+str(nr_variables))
    print('number of hospitals: '+str(maxSplits))
    print('save eigenvalues: '+str(save_eigen))
    print('repeats: '+str(repeat))
    print('added noise: '+str(noise))

    data, sample_ids, variable_names = data_import(datafile, ddpca, seed, nr_variables, header, rownames,  outfile= dirname)

    try:
        os.mkdir(dirname)
    except OSError:
        print("Creation of the directory %s failed" % dirname)
    if noise:
        if epsilons is None:
            epsilons = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 100]
        if deltas is None:
            deltas = [0.01]
    # projections
        filename = dirname + '/'
        res, val, vec = simul.run_multiple_simulations(data, dims=ndims, dppca=ddpca, noise=noise, nrSamples=repeat, splits=maxSplits,epsilons=epsilons,deltas=deltas, dirname=dirname, save_eigen=save_eigen, filename =filename)

    else:
        filename = dirname + '/nonoise_'
        res, val, vec = simul.run_multiple_simulations(data, dims=ndims, dppca=ddpca, noise=noise, nrSamples=repeat,
                                             splits=maxSplits, dirname=dirname,
                                             save_eigen=save_eigen, filename=filename)



def projection(scaled, sim, filename):
    projection = sc.dot(scaled, sim[:, 0:2])
    projection = pd.DataFrame(projection)
    projection.to_csv(filename, sep='\t', header=None, index=False)
#3def run_premade_datasets(datasets):
 #   for d in datasets:
  #      data, sample_ids, variable_names = data_import(d, ddpca, seed, nr_variables, header, rownames,outfile=dirname)


if __name__ == "__main__":

    parser = ap.ArgumentParser(description='Run distributed PCA simulation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='output file', type=str, help='filename of standalone PCA output')
    parser.add_argument('-d', metavar='dimensions', type=int, help='number of principal components to return')
    parser.add_argument('-p', metavar='output directory', type=str, help='output directory for simulation study.')
    parser.add_argument('-s', metavar='random seed', type=int, help='Seed for random processes')
    parser.add_argument('-n', metavar='Number of samples', type=int, help='Number of variables to sample in the case of high dimensional datasets', default=None)
    parser.add_argument('-k', metavar='Number hospitals',type=int, help='Number of simulated hospitals', default=5)
    parser.add_argument('-e', metavar='Save eigenvalues', type=bool, help='If true the generated eigenspaces are saved (!a lot of large files)', default=False)
    parser.add_argument('-r', metavar='Repeats', type=int, help='Number of repetitions of the sampling process', default=5)
    parser.add_argument('-c', metavar='Colnames', type=bool, help='True of data has column headers',
                        default=False)
    parser.add_argument('-i', metavar='Sample ids', type=int, help='Dataframe column which contains the sample ids, 0 index,',default=-1)
    args = parser.parse_args()


    print('file: '+str(args.f))
    print('output file: '+str(args.o))
    print('dimensions: '+str(args.d))
    print('output directory: '+str(args.p))
    print('random seed: '+str(args.s))
    print('number of sampled variables: '+str(args.n))
    print('number of hospitals: '+str(args.k))
    print('save eigenvalues: '+str(args.e))
    print('repeats: '+str(args.r))
    print('column headers: '+ str(args.c))
    print('sample ids: '+ str(args.i))

    ddpca = Distributed_DP_PCA()

    if args.c :
        header = 0
    else:
        header = None

    if args.i == -1:
        rownames=None
    else:
        rownames = args.i

    '''
    import scipy.linalg as la
    import scipy.sparse.linalg as lsa
    datafile ='/home/anne/Documents/featurecloud/data/TCGA/htseq/out.trunc.txt'
    out = '/home/anne/Documents/featurecloud/data/TCGA/htseq/out.text'
    outdir ='/home/anne/Documents/featurecloud/data/TCGA/htseq/'
'''
    ddpca = Distributed_DP_PCA()
    """
    data, d,s = data_import(datafile, ddpca, seed=11, nr_samples=None, header=None, rownames=4)
    u,v,w = la.svd(dat)
    w = pd.DataFrame(sc.transpose(w))
    w.to_csv(path_or_buf='/home/anne/Documents/featurecloud/data/TCGA/htseq/svd.eigen', sep='\t')
    dat = pd.DataFrame(dat)
    dat.to_csv(path_or_buf='/home/anne/Documents/featurecloud/data/TCGA/htseq/scaled.txt', sep='\t')
    run_standalone(data, ddpca, outfile=out, dims=2, seed=22, nr_samples=None)
	"""
    '''
    run_standalone(datafile, ddpca, outfile= out, dims=4, seed=11, header=0, rownames=0)
    run_1_n_sites_with_noise(datafile, ddpca, outdir, False, ndims=4, maxSplits=5, repeat=1, save_eigen = True, seed = 11, header=0, rownames=0)
    run_1_n_sites_with_noise(datafile, ddpca, outdir, True, ndims=4, maxSplits=5, repeat=5, save_eigen=True, seed=11, header=0, rownames=0)

    '''
    run_standalone(args.f, ddpca, outfile= args.o, dims=args.d, seed=args.s, nr_samples = args.n, header = header, rownames = rownames)
    run_1_n_sites_with_noise(args.f, ddpca, args.p, False, ndims=args.d, maxSplits=args.k, repeat=1,  save_eigen=args.e, seed=args.s, nr_variables=args.n,header = header, rownames = rownames)
    run_1_n_sites_with_noise(args.f, ddpca, args.p, True, ndims=args.d, maxSplits=args.k, repeat=args.r, save_eigen=args.e, seed=args.s, nr_variables=args.n,header = header, rownames = rownames)

