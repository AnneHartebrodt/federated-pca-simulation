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


def leave1out_site(data_list, outdir_name, dataset_name, k=20):
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')

    u, s, v = compute_canonical(data_list, k=k)
    data_list_copy = data_list[:]
    for i in range(len(data_list)):
        data_list.pop(i)
        u1, s1, v1 = compute_canonical(data_list, k=k)
        angles2 = co.compute_angles(v, v1, reported_angles=k)
        with open(filename, 'a') as handle:
            handle.write(cv.collapse_array_to_string(angles2, str(i)))
        data_list = data_list_copy[:]

def leave1out_random_subset(data_list, outdir_name, dataset_name, k=20):
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name + '_angles.tsv')

    # get all dataset sizes
    sizes = []
    for d in data_list:
        sizes.append(d.shape[0])

    data = np.concatenate(data_list, axis=0)
    u, s, v = compute_canonical(data, k=k)
    for i in range(10):
        for s in sizes:
            # generate a random sample of the size of one of the sites
            drop = np.random.choice(range(data.shape[0]), s, replace=False)
            u1, s1, v1 = compute_canonical(np.delete(data, drop, axis=0), k=k)
            angles2 = co.compute_angles(v, v1, reported_angles=k)
            with open(filename, 'a') as handle:
                handle.write(cv.collapse_array_to_string(angles2, str(s)+'\t'+str(i)))

def leave1out_random_sample(data_list, outdir_name, dataset_name, row_list, k=20):
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')

    data = np.concatenate(data_list, axis=0)
    row_list = np.concatenate(row_list)
    u, s, v = compute_canonical(data, k=k)

    # generate a random sample of 100 points to drop
    #drop = np.random.choice(range(data.shape[0]), 100, replace=False)
    drop = range(data.shape[0])
    for d in drop:
        u1, s1, v1 = compute_canonical(np.delete(data, d, axis=0), k=k)
        angles2 = co.compute_angles(v, v1, reported_angles=k)
        with open(filename, 'a') as handle:
            handle.write(cv.collapse_array_to_string(angles2, str(row_list[d])))

def leave1out_random_sample_balcan(data_list, outdir_name, dataset_name, row_list, k=20):
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')

    data = np.concatenate(data_list, axis=0)
    u, s, v = compute_canonical(data, k=k)

    # generate a random sample of 100 points to drop
    #drop = np.random.choice(range(data.shape[0]), 100, replace=False)
    # iterate over all lists
    for l in range(len(data_list)):
        # iterate over all samples in the list
        for i in range(data_list[l].shape[0]):
            cp = data_list[l].copy()
            data_list[l] = np.delete(data_list[l], i, axis=0)
            v1, s1= b.simulate_federated_horizontal_pca(data_list, k=k)
            data_list[l] = cp
            print(row_list[l][i])
            angles2 = co.compute_angles(v, v1, reported_angles=k)
            with open(filename, 'a') as handle:
                handle.write(cv.collapse_array_to_string(angles2, str(row_list[l][i])))

if __name__ == '__main__':
    basedir = '/home/anne/Documents/featurecloud/data/tcga/cancer_type/'
    datasets = os.listdir(basedir)
    data_dirs = [op.join(basedir, d, 'sites', 'data') for d in datasets]

    #datasets = [datasets[1]]
    #data_dirs = [data_dirs[1]]
    ## presplit test
    outdir_leave1site = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/site'
    outdir_leave1sample = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/sample'
    outdir_leave1site_random = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/random'
    outdir_leave1site_balcan = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/sample_balcan'
    k = 20
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            filenames = os.listdir(d)
            data_list, row_list = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            #u, s, v = compute_canonical(data_list, k)
            leave1out_site(data_list=data_list, outdir_name=outdir_leave1site, k=k, dataset_name=name)

            data_list, row_list = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            leave1out_random_sample(data_list, outdir_leave1sample,dataset_name=name, row_list=row_list, k=k)

            data_list, row_list = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            leave1out_random_subset(data_list=data_list, outdir_name=outdir_leave1site_random, k=k, dataset_name=name)

            data_list, row_list = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            leave1out_random_sample_balcan(data_list, outdir_leave1site_balcan ,dataset_name=name, row_list=row_list, k=k)
        except FileNotFoundError:
            print('File not found')


