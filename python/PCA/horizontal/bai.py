import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *


def combine_subroutine(data_list, m_list):
    """ avoid as many QR factorisations as possible to
    reduce numerical issues."""
    i = 0
    qr_list = []
    for d, m in zip(data_list, m_list):
        print('hello')
        q, r = la.qr(d, mode='economic')
        m = np.atleast_2d(np.nanmean(data, axis=0))
        #r = np.concatenate([r, m], axis=0)
        qr_list.append(r)
    g, gl, r, rl = simulate_federated_qr(qr_list)
    return r

from python.PCA.vertical.simulate_federated_qr_orthonormalisation import simulate_federated_qr
def simulate_bai(data_list, k=10):

    local_sums = [np.sum(d, axis=0) for d in data_list]
    sample_count = [d.shape[0] for d in data_list]
    local_means = [local_sums[i]/sample_count[i] for i in range(len(local_sums))]

    global_means = np.sum(local_sums, axis=0)
    global_means = global_means/(sum(sample_count))

    m = [np.sqrt(sample_count[i])* (local_means[i]-global_means[i]) for i in range(len(sample_count))]

    r= combine_subroutine(data_list, m)
    u, s, v = lsa.svds(r, k=k)
    v = np.flip(v.T, axis=1)
    s = np.flip(s)
    u = np.flip(u)
    ul = [np.dot(d,v) for d in data_list]
    return ul, s, v


if __name__ == '__main__':
    from python.PCA.horizontal.horizontal_pca_benchmark import read_presplit_data_folders, compute_canonical,scale_datasets


    # MNIST for reference
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    data_list, choices = sh.partition_data_horizontally(data, splits=4, randomize=False)
    u, s, v = simulate_bai(data_list, k=10)
    data = np.concatenate(data_list, axis=0)
    uu, ss, vv = lsa.svds(data, k=10)
    vv = np.flip(vv.T, axis=1)
    ang = co.compute_angles(vv, v)
