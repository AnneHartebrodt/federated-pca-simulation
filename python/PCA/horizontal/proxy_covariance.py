from python.PCA.vertical.vertical_pca_benchmark import *
import numpy as np
import os

def cross_site_variance(local_sums, sample_count):
    local_sums = np.float32(local_sums)
    local_means = [local_sums[i]/sample_count[i] for i in range(len(sample_count))]
    totaln = np.sum(sample_count)
    sum = np.sum(local_sums, axis=0)
    global_mean = sum/totaln
    # local means - global means > cross site covariance
    mm = np.atleast_2d(local_means[0] - global_mean)
    matr = sample_count[0]*np.dot(mm.T, mm)
    for m in range(1,len(local_sums)):
        mm = np.atleast_2d(local_means[m] - global_mean)
        matr = matr + sample_count[m] * np.dot(mm.T, mm)
    matr = matr/totaln
    return matr

def simulate_federated_horizontal_pca_Qu(datasets, local_sums, sample_count, k=10, qu=True, weighted=False, k_factor=2):
    partial = []
    for d in datasets:
        # Convert to float 32 to save memory
        d = np.float32(d)
        partial.append(local_SVD_Qu(d, k=k_factor*k))
    dpca = aggregate_partial_SVDs_Qu(partial, local_sums=local_sums, sample_count=sample_count, t1=k, qu=qu, weighted=weighted)
    return dpca

def aggregate_partial_SVDs_Qu(svd_list, local_sums, sample_count, t1=10, qu=True, weighted=False):
    """
    This function aggregates the local proxy covariances by averaging them

    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:

    """
    k1 = min(sample_count[0]-1, svd_list[0][0].shape[1])
    R = np.zeros((k1, k1))
    np.fill_diagonal(R, svd_list[0][1])
    total_count = np.sum(sample_count)
    if qu or weighted:
        Ac = sample_count[0] * np.dot(np.dot(svd_list[0][0], R), svd_list[0][0].T)
    else:
        Ac = np.dot(np.dot(svd_list[0][0], R), svd_list[0][0].T)
    for svd in range(1, len(svd_list)):
        k1 = min(sample_count[svd]-1, svd_list[svd][0].shape[1])
        R = np.zeros((k1, k1))
        np.fill_diagonal(R, svd_list[svd][1])
        if qu or weighted:
            Ac = Ac + sample_count[svd] * np.dot(np.dot(svd_list[svd][0], R), svd_list[svd][0].T)
        else:
            Ac = Ac + np.dot(np.dot(svd_list[svd][0], R), svd_list[svd][0].T)
    if qu:
        print('Adding cross site variance')
        Ac = Ac + cross_site_variance(local_sums, sample_count)

    if weighted:
        Ac = Ac/total_count

    U, S, UT, nd = sh.svd_sub(Ac, t1)
    S = np.sqrt(S)
    return UT, S


def local_SVD_Qu(data, k=20):
    """
    Performs a singular value decomposition local data
    :param cov: A covariance matrix
    :param r: The number of top principal components to be considered
    :return: U_r*S_r (The product of the matrices taking the top r colums/rows)
    """

    # returns column vectors only
    U, S, UT, nd = sh.svd_sub(data, k)
    # In case we want to use more values in the approximation

    return UT, S

def simulate_proxy_naive(data_list, k=10):
    samples = [d.shape for d in data_list]
    sample_count = np.sum(samples)
    proxy_covariance_matrix = np.dot(data_list[0].T, data_list[0])
    print(proxy_covariance_matrix.shape)
    for d in range(1, len(data_list)):
        proxy_covariance_matrix = proxy_covariance_matrix + np.dot(data_list[d].T, data_list[d])
    proxy_covariance_matrix = proxy_covariance_matrix/(sample_count-1)
    u,s,v, nd = sh.svd_sub(proxy_covariance_matrix, k)
    return u, s, v

if __name__ == '__main__':

    # # this is here to avoid circular import
    # from python.PCA.horizontal.horizontal_pca_benchmark import compute_canonical, scale_datasets, read_presplit_data_folders
    #
    #
    #
    # filenames = os.listdir('/home/anne/Documents/featurecloud/data/tcga/cancer_type/Liver_and_intrahepatic_bile_ducts/sites/data')
    # data_list, row_names = read_presplit_data_folders(filenames, '/home/anne/Documents/featurecloud/data/tcga/cancer_type/Liver_and_intrahepatic_bile_ducts/sites/data')
    # local_sums = [np.sum(d, axis=0) for d in data_list]
    # sample_count = [d.shape[0] for d in data_list]
    # # this algorithm works on the locally scaled data!!!!
    # # !!!
    # for d in range(len(data_list)):
    #     data_list[d] = scale_datasets([data_list[d]])[0]
    # x3, e = simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=True)
    #
    # dl = scale_datasets(data_list)
    # un, sn, vn = compute_canonical(dl, k=10)
    # ang = co.compute_angles(x3, vn)


    # MNIST for reference
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data_list, rows = sh.partition_data_horizontally(data, splits=3)


    # naive covariance method
    data_list = scale_datasets(data_list)
    un, sn, vn = simulate_proxy_naive(data_list, k=10)
    u, s, mv = compute_canonical(data_list, k=10)
    print(co.compute_angles(mv, vn))


    # resplit scaled data
    data_list, rows = sh.partition_data_horizontally(data, splits=3)
    # naive covariance method
    data_list = scale_datasets(data_list)
    local_sums = [np.sum(d, axis=0) for d in data_list]
    sample_count = [d.shape[0] for d in data_list]
    x4, e = simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=False)
    print(co.compute_angles(x4, mv))
    x4, e = simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=False, weighted=True)
    print(co.compute_angles(x4, mv))

    # MNIST for reference
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data_list, choice = sh.partition_data_horizontally(data, splits=2)

    # filenames = os.listdir(
    #    '/home/anne/Documents/featurecloud/data/tcga/cancer_type/Liver_and_intrahepatic_bile_ducts/sites/data')
    # data_list, row_names = read_presplit_data_folders(filenames,
    #                                                  '/home/anne/Documents/featurecloud/data/tcga/cancer_type/Liver_and_intrahepatic_bile_ducts/sites/data')

    local_sums = [np.sum(d, axis=0) for d in data_list]
    sample_count = [d.shape[0] for d in data_list]
    for d in range(len(data_list)):
        data_list[d] = scale_datasets([data_list[d]])[0]
    x3, e = simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=True)

    #ds = si.scale_center_data(data, center=True, scale_variance=False)
    #ds= scale_datasets(data_list)
    ds = scale_datasets(data_list)
    u, s, mv = compute_canonical(ds, k=10)
    print(co.compute_angles(x3, mv))






