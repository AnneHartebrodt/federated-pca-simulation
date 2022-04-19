import pandas as pd
import os
import scipy.sparse.linalg as lsa
from python.PCA.logging import *
import time
import os.path as op
import python.PCA.horizontal.horizontal_pca_power_iteration as h
import python.PCA.horizontal.balcan as b
import python.PCA.vertical.simulate_federated_vertically_partionned_pca as vertical
import python.import_export.mnist_import as mi
import python.PCA.shared_functions as sh
from scipy.sparse.coo import coo_matrix
import python.import_export.spreadsheet_import as si
from python.PCA.horizontal.horizontal_pca_benchmark import compute_canonical, read_presplit_data_folders, scale_datasets

def approximate_vertical(data_list, k=10, factor_k=2):
    v, e = b.simulate_federated_horizontal_pca(data_list, k=k, factor_k=factor_k)
    g = [np.dot(d, v) for d in data_list]
    return g


if __name__ == '__main__':
    # check if it performs on mnist
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    k = 20
    data_list, choices = sh.partition_data_horizontally(data, splits=20, randomize=False)
    data_list = scale_datasets(data_list)
    u, s, v = compute_canonical(data_list, k=20)

    g = approximate_vertical(data_list, 10, 2)

    g2 = np.concatenate(g, axis =0)

    co.compute_angles(u, g2)