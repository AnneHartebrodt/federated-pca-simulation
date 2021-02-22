
import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *



def simulate_federated_horizontal_pca(datasets, k=10):
    partial = []
    for d in datasets:
        partial.append(local_SVD(d, k=2*k))
    dpca = aggregate_partial_SVDs(partial, k)
    return dpca


def local_SVD(data, k=100):
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
    return P

def aggregate_partial_SVDs(svds, t2=100):
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

def simulate_federated_horizontal_pca_pc(datasets, k=10, means=None, weights=None, totaln=None):
    partial = []
    for d in datasets:
        partial.append(local_SVD(d, k=2*k))
    dpca = aggregate_partial_SVDs_pc(partial, k, means=means, weights=weights, totaln=totaln)
    return dpca

def cross_site_variance(means, weights, totaln):
    sum = np.sum(means, axis=0)
    sum = sum/totaln

    mm = np.reshape(means[0]/(totaln*weights[0]) - sum, (len(means[0]),1))
    matr =  np.dot(mm, mm.T) #*weights[0]
    for m in range(len(means)):
        mm = np.reshape(means[m]/(totaln*weights[m]) - sum, (len(means[m]),1))
        #matr = matr + (weights[m]* np.dot(mm, mm.T))
        matr = matr + (np.dot(mm, mm.T))
    return matr


def aggregate_partial_SVDs_pc(svd_list, t1=10, weights=None, means =None, totaln=None):
    """
    This function aggregates the local proxy covariances by averaging them

    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:

    """
    # Average covariance matrices
    s = len(svd_list)
    # by default we take all dimensions available
    if weights is not None:
        Ac = weights[0] * (np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :]))
        for svd in range(1, len(svd_list)):
            Ac = Ac + weights[svd] * (
                np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :]))
    else:
        Ac = np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :])
        for svd in range(1, len(svd_list)):
            Ac = Ac + np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :])
        Ac = 1 / s * Ac

    if means is not None:
        Ac  = Ac+ cross_site_variance(means, weights, totaln)

    U, S, UT, nd = sh.svd_sub(Ac, t1)
    S = np.sqrt(S)
    return UT[:, 0:nd], S[0:nd]


def read_presplit_data_folders(file_list, basedir):
    data_list = []
    for file in file_list:
        data = pd.read_csv(os.path.join(basedir,file), sep='\t')
        data = data.values
        data_list.append(data)
    return data_list

def compute_canonical(data_list, k=10):
    data = np.concatenate(data_list, axis=0)
    u,s,v = lsa.svds(data, k=k)
    u = np.flip(u, axis=1)
    v = np.flip(v.T, axis=1)
    s = np.flip(s)
    return u,s,v

def compute_equal_horizontal_split(data_list):
    data = np.concatenate(data_list, axis=0)
    data_list, choices = sh.partition_data_horizontally(data, splits=2, randomize=False)
    x, e = simulate_federated_horizontal_pca(data_list)
    return x,e

if __name__ == '__main__':

    basedir = '/home/anne/Documents/featurecloud/data/tcga/cancer_type/Breast/sites/data'
    filenames = os.listdir(basedir)
    data_list = read_presplit_data_folders(filenames,basedir)
    u,s,v = compute_canonical(data_list, k=10)
    x, e = simulate_federated_horizontal_pca(data_list)
    co.compute_angles(x, v)
    xx, ee = compute_equal_horizontal_split(data_list)
    co.compute_angles(xx,v)

    x2, e = simulate_federated_horizontal_pca_pc(data_list)
    co.compute_angles(x2, v)
    means = [np.sum(d, axis=0) for d in data_list]
    weights=[d.shape[0] for d in data_list]
    totaln= sum(weights)
    weights = [w/totaln for w in weights]
    x3, e = simulate_federated_horizontal_pca_pc(data_list, means=means, weights=weights, totaln=totaln)
    co.compute_angles(x3, v)




    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    #data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    data_list, choices = sh.partition_data_horizontally(data, splits=2, randomize=False)
    mx1, e = simulate_federated_horizontal_pca(data_list)
    u,s,mv = lsa.svds(data, k=10)
    co.compute_angles(mx1,np.flip(mv.T, axis=1))

    mx2, e = simulate_federated_horizontal_pca_pc(data_list)
    co.compute_angles(mx2, np.flip(mv.T, axis=1))

    means = [np.sum(d, axis=0) for d in data_list]
    weights=[data.shape[0] for d in data_list]
    totaln= 600000
    mx3, e = simulate_federated_horizontal_pca_pc(data_list, means=means, weights=weights, totaln=totaln)
    co.compute_angles(mx3, np.flip(mv.T, axis=1))
