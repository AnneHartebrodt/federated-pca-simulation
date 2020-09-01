import numpy as np
import guo_vertical as gv
import shared_functions as sh
import scipy.linalg as la
import import_export.easy_import as easy
import import_export.import_data as imp
import scipy.sparse.linalg as lsa
import comparison as co
import argparse as ap
import pandas as pd
import os.path as path



def simulate_guo(local_data, k, maxit):
    G_list = []
    iterations = 0
    ra = False
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
    return G_i

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='infile', type=str, help='filename of data file; default tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='filename of data file; default tab separated')
    parser.add_argument('-g', metavar='grm', type=str, default=None)
    parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
    parser.add_argument('-s', metavar='sites', default=10, type=int, help='Number of sites simulated')
    parser.add_argument('-i', metavar='iteration', default=2000, type=int, help='Maximum number of iterations')
    parser.add_argument('-p', metavar='outpath', type=str, help='Output directory for result files')
    args = parser.parse_args()

    args.f =  '/home/anne/Documents/featurecloud/gwas/data/hapmap/thin.rec.raw.T.scaled.man'
    args.p = '/home/anne/Documents/featurecloud/gwas/results/hapmap'
    # import scaled SNP file
    data = easy.easy_import(args.f, header=None, rownames=None, center=False, scale_var=False,sep='\t')

    g = gv.standalone(data, k=args.k)

    u, s, v = lsa.svds(data.T,k=args.k)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)


    data_list = sh.partition_data_vertically(data,args.s)
    ug = simulate_guo(data_list,args.k+2, maxit=2000)


    pd.DataFrame(g).to_csv(path.join(args.p,'guo_single_site_eigenvector.tsv'), sep = '\t', header=False, index=False)
    pd.DataFrame(u).to_csv(path.join(args.p,'scipy_single_site_eigenvector.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(ug).to_csv(path.join(args.p,'guo_multi_site_eigenvector.tsv'), sep='\t', header = False, index=False)




