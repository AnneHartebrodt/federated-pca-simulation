import numpy as np
import guo_vertical as gv
import shared_functions as sh
import scipy.linalg as la
import import_export.easy_import as easy
import import_export.import_data as imp
import scipy.sparse.linalg as lsa
import comparison as co



def simulate_guo(local_data, k, maxit):
    G_list = []
    iterations = 0
    ra = False
    #dim = min(local_data[0].shape[1], k)
    #for d in local_data:
    #    dim = np.min(min(d.shape[1], dim))
    #for d in local_data:
    #    G_i = sh.generate_random_gaussian(d.shape[1], dim)  # phi1
    #    G_i, Q = la.qr(G_i, mode='economic')
    #    G_list.append(G_i)
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    while not ra and iterations<maxit:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], k))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            #G_i, Q = la.qr(G_i, mode='economic')
            G_list_n.append(G_i)
        start = 0
        G_i = np.concatenate(G_list_n, axis=0)
        G_i, R = la.qr(G_i, mode='economic')
        G_list = []
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]
        ra = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
    G_i = np.concatenate(G_list)
    return G_i, G_list

def simulate_guo_2(local_data):
    G_list = []
    iterations = 0
    ra = False
    for d in local_data:
        G_i = sh.generate_random_gaussian(d.shape[1], 1)  # phi1
        G_i = G_i / np.linalg.norm(G_i)
        G_list.append(G_i)
    H_i_prev = np.ones((local_data[0].shape[0],1))
    while not ra:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0],1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            G_i = G_i/np.linalg.norm(G_i)
            G_list_n.append(G_i)
        ra = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
        G_list = G_list_n
    G_i = np.concatenate(G_list)
    return G_i, G_list


if __name__ == '__main__':

    file = '/home/anne/Documents/data/mnist/flat.csv'
    outfile = '/home/anne/Documents/featurecloud/gwas/chr10'
    header = 0
    rownames = None
    center = False
    scale = False
    scale_var = False
    scale01 = False
    scale_unit = False
    p = 23
    transpose = False
    sep = ','
    #
    import pandas as pd

    data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                            scale01=scale01, scale_unit=scale_unit,
                            outfile=outfile, transpose=transpose, sep=sep)
    #data = imp.CustomDataImporter().scale_center_data(data, center=True)
    #data = data[:,0:(data.shape[1]-1)]

    g = gv.standalone(data, k=2)

    u, s, v = lsa.svds(data.T,k=2)
    u = np.flip(u, axis=1)

    co.compute_angles(g, u)

    data_list = sh.partition_data_vertically(data,2)

    r, l = simulate_guo_2(data_list)
    co.angle(r.flatten(), u[:,0])

    r, l = simulate_guo(data_list,1)
    co.compute_angles(r, u)


