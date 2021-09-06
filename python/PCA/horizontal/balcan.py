import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *


def simulate_federated_horizontal_pca(datasets, k=10, factor_k=2):
    partial = []
    for d in datasets:
        #d = np.float32(d)
        partial.append(local_SVD(d, k=factor_k*k))
    print('Intermediate dimensions' +str(factor_k*k))
    dpca = aggregate_partial_SVDs(partial, factor_k*k)
    return dpca


def local_SVD(data, k=20):
    """
    Performs a singular value decomposition local data
    :param cov: A covariance matrix
    :param r: The number of top principal components to be considered
    :return: U_r*S_r (The product of the matrices taking the top r colums/rows)
    """

    # returns column vectors only
    U, S, UT, nd = sh.svd_sub(data, k)
    # In case we want to use more values in the approximation
    nd = min(nd, k)
    R = np.zeros((nd, nd))
    np.fill_diagonal(R, S[0:nd])
    U_r = UT[:, 0:nd]
    P = np.dot(np.sqrt(R), U_r.T)
    print(P.shape)
    return P

def aggregate_partial_SVDs(svds, t2=10):
    """
    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:
    """

    svds = np.concatenate(svds, axis=0)
    ndim = min(t2, svds.shape[0] - 1)

    print(svds.shape)
    U, S, UT, nd = sh.svd_sub(svds, ndim)
    return UT, S



if __name__ == '__main__':
    # this is here to avoid circular import
    from python.PCA.horizontal.horizontal_pca_benchmark import read_presplit_data_folders, compute_canonical, scale_datasets


    # MNIST for reference
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)

    #ds = si.scale_center_data_columnwise(data, center=True, scale_variance=True)

    data_list, choice = sh.partition_data_horizontally(data, splits=2, randomize=False)

    data_list = scale_datasets(data_list)
    sample_count = [d.shape[0] for d in data_list]
    total_samples = sum(sample_count)
    weights = [sc/total_samples for sc in sample_count]
    u,s,v = compute_canonical(data_list, k=20)
    uu, ss, vv = lsa.svds(np.concatenate(data_list, axis=0))
    x, e = simulate_federated_horizontal_pca(data_list, k=716, factor_k=3)
    x2, e = simulate_federated_horizontal_pca(data_list, k=20, factor_k=2)
    v1 = v
    v1[0,0]=-8.36060474094277e-05
    v1[5,0] =-0.00056577135205741441
    co.compute_angles(x, v)[12]
    np.sum(x[:,12]-v[:,12])/717
    co.compute_angles(x, v)
    print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), v))
    #print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), vv))
    print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), x))


