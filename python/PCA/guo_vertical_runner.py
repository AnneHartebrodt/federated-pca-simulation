'''
    Copyright (C) 2020 Anne Hartebrodt

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    Authors: Anne Hartebrodt

'''

import numpy as np
import python.PCA.guo_vertical as gv
import python.PCA.shared_functions as sh
import scipy.linalg as la
#import import_export.easy_import as easy
import scipy.sparse.linalg as lsa
import argparse as ap
import pandas as pd
import os.path as path
import python.import_export.mnist_import as imnist
import python.PCA.comparison  as co
from scipy.sparse import coo_matrix
import os


def simulate_guo(local_data, k, maxit):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
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
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        G_i, R = la.qr(G_i, mode='economic')
        G_list = []
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]
        ra = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i

    G_i = np.concatenate(G_list)
    return G_i, eigenvals

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
    from scipy.sparse import coo_matrix
    # import scaled SNP file
    #data = easy.easy_import(args.f, header=None, rownames=None, center=False, scale_var=False,sep='\t')
    data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    data = coo_matrix.asfptype(data)
    args.k = 10
    g = gv.standalone(data, k=10)

    u, s, v = lsa.svds(data.T,k=10)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)


    data_list, choices = sh.partition_data_vertically(data,2)
    ug, ev = simulate_guo(data_list, 12, maxit=200)

    co.compute_angles(ug, u)

    pd.DataFrame(g).to_csv(path.join(args.p,'guo_single_site_eigenvector.tsv'), sep = '\t', header=False, index=False)
    pd.DataFrame(u).to_csv(path.join(args.p,'scipy_single_site_eigenvector.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(ug).to_csv(path.join(args.p,'guo_multi_site_eigenvector.tsv'), sep='\t', header = False, index=False)




