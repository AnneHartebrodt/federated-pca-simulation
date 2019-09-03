from PCA.master import Distributed_DP_PCA
import pandas as pd
import scipy as sc
import scipy.spatial.distance as d
import numpy as np
import random as r
import copy
import os

def calculate_euclidean_distance(V1, V2):
    res = []
    for line in range(min(V1.shape[1], V2.shape[1])):
        res.append(d.euclidean(V1[:,line], V2[:,line]))
    return (res)

def run_multiple_simulations(dataset, dims, noise, splits, nrSamples, dppca, epsilons=[0.01], deltas=[0.01], dirname='eigenspaces', save_eigen=False, filename = ''):
    # Just to be sure.
    print('Running simulatiom with noise: '+str(noise))
    backup = copy.deepcopy(dataset)
    dm = min(dims, dataset.shape[1])
    results = np.empty(shape=(1,dm+4))
    eigenvalues = np.empty(shape=(1,dm+4))
    eigenvectors = np.empty(shape=(1,dm+4))
    for epsilon in epsilons:
        for delta in deltas:
            v1 = None
            for i in range(nrSamples):
                for split in range(1, splits + 1):
                    if noise:
                        currentFile = dirname+'/eigenvectors_split-'+str(split)+'_i-'+str(i)+'_epsilon-'+str(epsilon)+'_delta-'+str(delta)+'.txt'
                    else:
                        currentFile = dirname + '/eigenvectors_split-' + str(split) + '_i-' + str(
                            i) + '_no_noise' + '.txt'
                    dataset = copy.deepcopy(backup)
                    vec,val = dppca.simulate_multisite_PCA(dataset, split, epsilon, delta, noise=noise, ndims=dims, scale=False, center=False)
                    if v1 is None:
                        v1 = copy.deepcopy(vec)
                    else:
                        vec = dppca.normalize_eigenspaces([v1,vec])[1]
                    #projection matrix
                    proj = sc.dot(dataset, vec[:, 0:dm])
                    proj = np.concatenate((proj, create_annotation(proj.shape[0], split, i, epsilon, delta)), axis=1)
                    results = np.concatenate((results, proj), axis=0)
                    # eigenvalues
                    ar = np.array(np.concatenate((val, np.array([split, i, epsilon, delta])))).reshape((1,dm+4))
                    eigenvalues = np.concatenate((eigenvalues, ar), axis=0)
                    #eigenvectors
                    vec = np.concatenate((vec, create_annotation(vec.shape[0], split, i, epsilon, delta)), axis=1)
                    eigenvectors = np.concatenate((eigenvectors, vec), axis=0)

    # remove first, random row
    results = np.delete(results, 0, 0)
    eigenvectors = np.delete(eigenvectors, 0,0)
    eigenvalues = np.delete(eigenvalues, 0, 0)

    pd.DataFrame(results).to_csv(filename+'projections.txt', header=False, index=False, sep='\t')
    pd.DataFrame(eigenvalues).to_csv(filename+'eigenvalues.tsv', header=False, index=False, sep='\t')
    pd.DataFrame(eigenvectors).to_csv(filename+'eigenvectors.tsv', header=False, index=False, sep='\t')

    return results, eigenvalues, eigenvectors


def create_annotation(n, split, i, epsilon, delta):
    proj = np.ones((n, 1)) * split
    proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * i)), axis=1)
    proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * epsilon)), axis=1)
    proj = np.concatenate((proj, (np.ones((proj.shape[0], 1)) * delta)), axis=1)
    return proj

def prepare_data_iris():
    original = pd.read_csv(
        filepath_or_buffer='https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data',
        header=None,
        sep=',')
    r.seed(a=12)
    original = original.sample(frac=1, random_state=12).reset_index(drop=True)  # randomize dataframe to get rid of the order
    original.columns = ['sepal_len', 'sepal_wid', 'petal_len', 'petal_wid', 'class']
    original.dropna(how="all", inplace=True)  # drops the empty line at file-end
    labels = original.iloc[:, 4].values
    orig = sc.array(original)[:, 0:4]
    orig = orig.astype(float)
    return orig, labels

def perform_standalone_pca(data, ndims, dppca, center=True, scale=True, scale01=False):
    """
    This function performs a standard principal component analysis via eigendecomposition
    of the covariance matrix
    :param data: input numpy data nxd with n the number of samples and d the number of variable
    :param ndims: Number of dimensions to project onto
    :param center: if true, data is centered
    :param scale: if true, data is scaled to equal variance
    :param scale01: if true, data is scaled between 0 and 1
    :return: the projected data, The first ndim eigenvectors and eigenvalues
    """
    data = dppca.scale_data(data, center, scale, scale01)
    n = data.shape[0] # get the number of rows
    cov = (1 / (n - 1)) * sc.dot(data.transpose(), data) # use unbiased estimator
    #print(cov)
    # the sparse matrix version of svd has better memory requirements, while being a little
    # slower
    # covariance matrix is positive semi definite so SVD= Eigenvalue decomposition
    nd = min(data.shape[1], ndims)
    #print(nd)
    if(nd>10):
        print('Using sparse PCA for decreased memory consumption')
        nd = min(data.shape[1] - 1, ndims)
        V_global, S, W = sc.sparse.linalg.svds(cov, nd)
    else:
        nd = min(data.shape[1], ndims)
        V_global, S, W = sc.linalg.svd(cov)

    W = np.transpose(W) # eigenvectors
    W = dppca.normalize_eigenvectors(W)
    # this changes nothing

    # create projection matrix by multiplying by the first ndim eigenvalues
    proj_global = sc.dot(data, W[:, 0:nd])
    return (proj_global, W[:, 0:ndims], S[0:nd])




if __name__=="__main__":
    print('main')

