from operator import index

import pandas as pd
import os
import scipy.sparse.linalg as lsa
from python.PCA.logging import *
import time
import os.path as op
import python.PCA.horizontal.horizontal_pca_power_iteration as h
import python.PCA.horizontal.balcan as b
import python.PCA.horizontal.bai as bai
import python.PCA.horizontal.proxy_covariance as proxy
import python.PCA.vertical.simulate_federated_vertically_partionned_pca as vertical
import python.import_export.mnist_import as mi
import python.PCA.shared_functions as sh
from scipy.sparse.coo import coo_matrix

import python.import_export.spreadsheet_import as si

def read_presplit_data_folders(file_list, basedir):
    data_list = []
    row_list = []
    for file in file_list:
        data = pd.read_csv(os.path.join(basedir,file), sep='\t', index_col=0)
        rows = data.index
        row_list.append(rows)
        data = data.values
        data_list.append(data)
    return data_list, row_list

def compute_canonical(data_list, k=10):
    if isinstance(data_list, list):
        data = np.concatenate(data_list, axis=0)
    else:
        data = data_list
    u,s,v = lsa.svds(data, k=k)
    u = np.flip(u, axis=1)
    v = np.flip(v.T, axis=1)
    s = np.flip(s)
    return u,s,v


def wrapper_k_variation(data_list, outdir_name, k,  min_factor_k=2, max_factor_k=5, dataset_name="simulation"):
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_sre = op.join(outdir, dataset_name+'_sre.tsv')


    for ik in range(min_factor_k, max_factor_k):
        u, s, precomputed_eigenvector = compute_canonical(data_list, k)

        # balcan version
        #print ('intermediate dimenions: '+ str(ik))
        xx, ee = b.simulate_federated_horizontal_pca(data_list, k=k, factor_k=ik)
        angles2 = co.compute_angles(xx, precomputed_eigenvector, reported_angles=k)
        sre2 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), xx)
        with open(filename, 'a') as handle:
            handle.write(cv.collapse_array_to_string(angles2, 'balcan\t'+str(ik)))
        with open(filename_sre, 'a') as handle:
            handle.write(cv.collapse_array_to_string(sre2, 'balcan\t'+str(ik)))

        # proxy covaraince
        sample_count = [d.shape[0] for d in data_list]
        local_sums = [np.sum(d, axis=0) for d in data_list]
        xx, ee = proxy.simulate_federated_horizontal_pca_Qu(data_list, local_sums, sample_count,
                                                            k=k, k_factor=ik, weighted=False, qu=False)
        angles2 = co.compute_angles(xx, precomputed_eigenvector, reported_angles=k)
        sre2 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), xx)

        with open(filename, 'a') as handle:
            handle.write(cv.collapse_array_to_string(angles2, 'proxy\t'+str(ik)))
        with open(filename_sre, 'a') as handle:
            handle.write(cv.collapse_array_to_string(sre2, 'proxy\t'+str(ik)))



def wrapper(data_list, outdir_name, precomputed_eigenvector, k=20, dataset_name="simulation", maxit=2000):
    start = time.monotonic()
    # run power iteration benchmark
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')
    filename_sre = op.join(outdir, dataset_name+'_sre.tsv')


    # standard distributed power iteration
    x, e, count = h.simulate_distributed_horizontal(data_list, k, maxit=maxit)
    angles = co.compute_angles(x, precomputed_eigenvector)
    sre = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), x)
    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles, 'power_iteration'))
    with open(filename_iteration, 'a') as handle:
        handle.write('power_iteration' +'\t' + str(count)+'\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre, 'power_iteration\t' + str(count)))

    # balcan version
    xx, ee = b.simulate_federated_horizontal_pca(data_list, k)
    angles2 = co.compute_angles(xx, precomputed_eigenvector)
    sre2 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), xx)

    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles2, 'balcan_proxy'))
    with open(filename_iteration, 'a') as handle:
        handle.write('balcan_proxy' +'\t' + str(1)+'\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre2, 'balcan_proxy\t' + str(1)))

    # bal
    # vertical power iteration
    data_list_T = [d.T for d in data_list]
    u, s, vt, count = vertical.simulate_subspace_iteration(data_list_T, k, maxit=maxit)
    angles2 = co.compute_angles(vt, precomputed_eigenvector)
    sre2 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), vt)

    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles2, 'vertical_pca'))
    with open(filename_iteration, 'a') as handle:
        handle.write('vertical_pca' +'\t' + str(count)+'\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre2, 'vertical_pca\t' + str(count)))

    # bal

    #bai version
    u1, s1, v1 = bai.simulate_bai(data_list, k=k)
    angles3 = co.compute_angles(v1, precomputed_eigenvector)
    sre3 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), v1)

    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(angles3, 'bai_qr'))
    with open(filename_iteration, 'a') as handle:
        handle.write('bai_qr' +'\t' + str(1)+'\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre3, 'bai_qr\t' + str(count)))


    local_sums = [np.sum(d, axis=0) for d in data_list]
    sample_count = [d.shape[0] for d in data_list]
    x4, e = proxy.simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=False, k=k, weighted=False)
    ang4 = co.compute_angles(x4, precomputed_eigenvector)
    sre4 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), x4)

    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(ang4, 'proxy'))
    with open(filename_iteration, 'a') as handle:
        handle.write('proxy' + '\t' + str(1) + '\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre4, 'proxy\t' + str(1)))

    # x4, e = proxy.simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=False, k=k, weighted=True)
    # ang2 = co.compute_angles(x4, precomputed_eigenvector)
    # with open(filename, 'a') as handle:
    #     handle.write(cv.collapse_array_to_string(ang2, 'proxy_weighted'))
    # with open(filename_iteration, 'a') as handle:
    #     handle.write('proxy_weighted' + '\t' + str(1) + '\n')

    # proxy covariance naive
    un, sn, vn = proxy.simulate_proxy_naive(data_list, k=k)
    ang3 = co.compute_angles(vn, precomputed_eigenvector)
    sre3 = co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), vn)

    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(ang3, 'proxy_naive'))
    with open(filename_iteration, 'a') as handle:
        handle.write('proxy_naive' + '\t' + str(1) + '\n')
    with open(filename_sre, 'a') as handle:
        handle.write(cv.collapse_array_to_string(sre3, 'proxy_naive\t' + str(1)))


def wrapper_qi(data_list, outdir_name, precomputed_eigenvector, k=20, dataset_name="simulation", maxit=2000):
    # run power iteration benchmark
    outdir = op.join(outdir_name, dataset_name)
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, dataset_name+'_angles.tsv')
    filename_iteration = op.join(outdir, dataset_name+'_iterations.tsv')

    local_sums = [np.sum(d, axis=0) for d in data_list]
    sample_count = [d.shape[0] for d in data_list]
    for d in range(len(data_list)):
        data_list[d] = scale_datasets([data_list[d]])[0]

    x4, e = proxy.simulate_federated_horizontal_pca_Qu(data_list, local_sums=local_sums, sample_count=sample_count, qu=True, k=k)
    ang2 = co.compute_angles(x4, precomputed_eigenvector)
    with open(filename, 'a') as handle:
        handle.write(cv.collapse_array_to_string(ang2, 'proxy_qi'))
    with open(filename_iteration, 'a') as handle:
        handle.write('proxy_qi' + '\t' + str(1) + '\n')

def merge_randomize_and_split(data_list, splits=2):
    data = np.concatenate(data_list)
    data_list, choices = sh.partition_data_horizontally(data, splits=splits, randomize=True, equal=True)
    return data_list, choices

def scale_datasets(data_list):
    sums = []
    sample_count = 0

    # mean
    sums = [np.sum(d, axis = 0) for d in data_list]
    sums = np.sum(sums, axis=0)
    sample_count = [d.shape[0] for d in data_list]
    total_count = sum(sample_count)
    means = [s/total_count for s in sums]
    for i in range(len(data_list)):
        data_list[i] = data_list[i] - means

    #variance

    vars = [np.sum(np.square(d), axis=0) for d in data_list]
    vars = np.sum(vars, axis = 0)
    vars = vars/(total_count-1)
    # variance  = 0
    delete = np.where(vars==0)
    vars = np.delete(vars, delete)
    for i in range(len(data_list)):
        data_list[i] = np.delete(data_list[i], delete, axis=1)

    for i in range(len(data_list)):
        data_list[i] = data_list[i]/np.sqrt(vars)
    return data_list

def greedy_bin(sizes, bins=3):
    '''
    Very simplistic greedy heuristic to pack the numbers into bins.
    Does not produce very good results, notably when the range of
    values to be binned is large, but sufficient for the time being.
    Returns: the bins and the original indices of the values in the bins.

    '''
    list_indices = [[] for x in range(bins)]
    sums = [0 for x in range(bins)]
    min_index = 0
    buckets = [[] for x in range(bins)]
    for s in range(len(sizes)):
        buckets[min_index].append(sizes[s])
        list_indices[min_index].append(s)
        sums[min_index] += sizes[s]
        min = np.inf
        for i in range(len(sums)):
            if sums[i] < min:
                min_index = i
                min = sums[i]
    return buckets, list_indices

def merge_sites(data_list, indices):
    sample_count = {}
    total = 0
    for d in range(len(data_list)):
        sample_count[d] = {'start':total , 'end':total+data_list[d].shape[0]}
        total = total + data_list[d].shape[0]

    new_data_list = []
    order = []
    for i in indices:
        ilist = []
        for k in i:
            ilist.append(data_list[k])
            order += list(range(sample_count[k]['start'], sample_count[k]['end']))
        ilist = np.concatenate(ilist, axis=0)
        new_data_list.append(ilist)
    return new_data_list, order

def presplit(datasets, data_dirs, outdir_presplit, k):
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            filenames = os.listdir(d)
            data_list, row_list = read_presplit_data_folders(filenames, d)


            data_list = scale_datasets(data_list)
            u, s, v = compute_canonical(data_list, k)
            wrapper(data_list, outdir_presplit, precomputed_eigenvector=v, k=k, dataset_name=name)

            # reload the data to have it unscaled
            # ugly I know. Deep copy better?
            #data_list, row_list = read_presplit_data_folders(filenames, d)
            # scaling inside the function.
            #wrapper_qi(data_list, outdir_presplit, precomputed_eigenvector=v, k=k, dataset_name=name)

        except FileNotFoundError:
            print('File not found')

def presplit_merge(datasets, data_dirs, outdir_presplit_merged, outdir_k_variation,k=20 ):
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            filenames = os.listdir(d)
            data_list, row_list  = read_presplit_data_folders(filenames, d)

            # determine sized of data sets and store
            sizes = []
            for d in data_list:
                sizes.append(d.shape[0])
            sizes = np.array(sizes)

            data_list = scale_datasets(data_list)
            u, s, v = compute_canonical(data_list, k)

            # Merge sites together
            for s in [2,5]:
                # if there are not enough datasets to reasonable make larger sites
                # skip.
                if len(data_list) < s*2:
                    continue
                bins, indices = greedy_bin(sizes, s)
                data_list_new, choices = merge_sites(data_list, indices)
                outdir_s = op.join(outdir_presplit_merged, str(s))
                wrapper(data_list_new, outdir_s, precomputed_eigenvector=v, k=k, dataset_name=name)

                # k variation
                outdir_k_variation_s = op.join(outdir_k_variation, str(s))
                wrapper_k_variation(data_list_new, outdir_name=outdir_k_variation_s, k=k, min_factor_k=2, max_factor_k=5, dataset_name=name)

            #choice = np.random.choice(len(sizes), len(sizes), replace=False)
        except FileNotFoundError:
            print('File not found')

def k_variation(datasets, data_dirs, outdir_k_var, k=10):
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            filenames = os.listdir(d)
            data_list, row_list  = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            wrapper_k_variation(data_list, outdir_k_var, k=k, min_factor_k=1, max_factor_k=5, dataset_name=name)
        except FileNotFoundError:
            print('File not found')

def save_scaled_data(datasets, data_dirs, outdir, basedir):
    for d, name in zip(data_dirs, datasets):
        try:
            # get all file names
            outdir_scaled = op.join(basedir, name, 'sites', 'data_all')
            os.makedirs(outdir_scaled, exist_ok=True)
            filenames = os.listdir(d)
            data_list, row_list = read_presplit_data_folders(filenames, d)
            data_list = scale_datasets(data_list)
            data = np.concatenate(data_list, axis=0)
            rows = np.concatenate(row_list)
            pd.DataFrame(data).to_csv(op.join(outdir_scaled, name+'.tsv'), header=False, index=False, sep='\t')
            sample_provenance = np.concatenate(row_list)
            pd.DataFrame(sample_provenance).to_csv(op.join(outdir_scaled, name + '.names.tsv'), header=False, index=False, sep='\t')
        except FileNotFoundError:
            print('File not found')

if __name__ == '__main__':


    outdir_presplit = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/pre_split'
    outdir = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/scaled_data'
    outdir_k_var = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/k_variation/pre_split'
    outdir_presplit_merged = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/accuracy/merged'
    outdir_k_variation_merged = '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/k_variation/merged'

    # check if it performs on mnist
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    k = 20
    data_list, choices = sh.partition_data_horizontally(data, splits=20, randomize=False)
    data_list = scale_datasets(data_list)
    u, s, v = compute_canonical(data_list, k=20)
    wrapper(data_list, outdir_presplit, precomputed_eigenvector=v, k=k, dataset_name='mnist')
    wrapper_k_variation(data_list, outdir_k_var, k=k, min_factor_k=1, max_factor_k=5, dataset_name='mnist')


    # the real deal
    basedir = '/home/anne/Documents/featurecloud/data/tcga/cancer_type/'
    datasets = os.listdir(basedir)
    data_dirs = [op.join(basedir, d, 'sites', 'data') for d in datasets]

    # save scaled data
    #save_scaled_data(datasets, data_dirs, outdir, basedir)

    ## presplit test
    k = 20
    presplit(datasets, data_dirs, outdir_presplit, k)
    k_variation(datasets, data_dirs, outdir_k_var, k)
    presplit_merge(datasets, data_dirs, outdir_presplit_merged, outdir_k_variation_merged, k)

