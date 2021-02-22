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




def wrapper(data_list, outdir_name, precomputed_eigenvector, k=20, dataset_name="simulation", maxit=2000):
    start = time.monotonic()
    # run power iteration benchmark
    # simultaneous
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')


    # standard distributed power iteration
    x, e, count = h.simulate_distributed_horizontal(data_list, k, maxit=maxit)
    angles = co.compute_angles(x, precomputed_eigenvector)
    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles, 'power_iteration'))
    with open(filename_iteration, 'a') as handle:
        handle.write('power_iteration' +'\t' + str(count)+'\n')

    # balcan version
    xx, ee = b.simulate_federated_horizontal_pca(data_list, k)
    angles2 = co.compute_angles(xx, precomputed_eigenvector)
    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles2, 'balcan_proxy'))
    with open(filename_iteration, 'a') as handle:
        handle.write('balcan_proxy' +'\t' + str(1)+'\n')


    data_list_T = [d.T for d in data_list]
    u, s, vt, count = vertical.simulate_subspace_iteration(data_list_T, k, maxit=maxit)
    angles2 = co.compute_angles(vt, precomputed_eigenvector)
    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles2, 'vertical_pca'))
    with open(filename_iteration, 'a') as handle:
        handle.write('vertical_pca' +'\t' + str(count)+'\n')

def merge_randomize_and_split(data_list, splits=2):
    data = np.concatenate(data_list)
    data_list, choices = sh.partition_data_horizontally(data, splits=splits, randomize=True, equal=True)
    return data_list, choices


if __name__ == '__main__':
    basedir = '/home/anne/Documents/featurecloud/data/tcga/cancer_type/'
    datasets = os.listdir(basedir)[1:2]
    data_dirs = [op.join(basedir, d, 'sites', 'data') for d in datasets]
    outdir_presplit = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/pre_split'
    k = 20
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            filenames = os.listdir(d)
            data_list = read_presplit_data_folders(filenames, d)
            u, s, v = compute_canonical(data_list)
            wrapper(data_list, outdir_presplit, precomputed_eigenvector=v, k=20, dataset_name=name)
        except FileNotFoundError:
            print('File not found')


    # do the same, but with radomized data
    for d, name in zip(data_dirs, datasets):
        # get all file names
        filenames = os.listdir(d)
        data_list = read_presplit_data_folders(filenames, d)
        # compute canonical pca
        u, s, v = compute_canonical(data_list)
        for splits in [2, 5]:
            outdir_randomized = op.join('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/randomized', str(splits))
            try:
                #randomize data
                data_list, choices = merge_randomize_and_split(data_list, splits=splits)
                # reorder precomputed eigenvecotr in the new order
                #precomputed = v[choices, :]
                wrapper(data_list, outdir_randomized, precomputed_eigenvector=v, k=20, dataset_name=name)
            except FileNotFoundError:
                print('File not found')

    # check if it performs on mnist
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    data_list, choices = sh.partition_data_horizontally(data, splits=20, randomize=False)
    u,s,v = compute_canonical(data_list, k)
    wrapper(data_list, outdir_presplit, precomputed_eigenvector=v, k=20, dataset_name='mnist')
