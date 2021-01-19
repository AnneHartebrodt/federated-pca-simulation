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
import python.PCA.vertical_pca_library as gv
import python.PCA.shared_functions as sh
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import python.import_export.mnist_import as imnist
import python.PCA.comparison  as co
from scipy.sparse import coo_matrix
import time

def simulate_guo_k(local_data, maxit, V_k):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.
    by imposing orthogonality via qr factorisation.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = [] # this are the parital eigevenctors
    iterations = 0
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    k = V_k.shape[1]
    G_i = sh.generate_random_gaussian(total_len, 1)
    G_i = np.concatenate([V_k, G_i], axis=1)
    G_i, R = la.qr(G_i, mode='economic')
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1],  k:])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not ra and iterations<maxit:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], 1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            G_list_n.append(G_i)
        G_i = np.concatenate(G_list_n, axis=0)
        G_i = np.concatenate([V_k, G_i], axis=1)
        G_i, R = la.qr(G_i, mode='economic')
        G_list = []
        start = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1],  k:])
            start = start+d.shape[1]
        ra, sum = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i

    G_i = np.concatenate(G_list)
    return G_i

def simulate_guo_1(local_data, maxit):
    '''
    Retrieve the first eigenvector
    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = [] # this are the parital eigevenctors
    iterations = 0
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, 1)
    G_i = G_i/np.linalg.norm(G_i) #normalize
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not ra and iterations<maxit:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], 1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            G_list_n.append(G_i)
        G_i = np.concatenate(G_list_n, axis=0)
        G_i = G_i/la.norm(G_i)
        G_list = []
        start = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start+d.shape[1]
        ra, sum = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i

    G_i = np.concatenate(G_list)
    return G_i


def simulate_guo_k_residuals_local(local_data, maxit, V_k):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = [] # this are the parital eigevenctors
    iterations = 0
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = residuals(V_k).T
    Vk_list = []
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        Vk_list.append(V_k[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not ra and iterations<maxit:
        iterations = iterations + 1

        H_i = np.zeros((local_data[0].shape[0], 1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        local_sums = []
        for d, g, v in zip(local_data, G_list, Vk_list):
            G_i = np.dot(d.T, H_i) + g
            G_list_n.append(G_i)
            sum = []
            for vi in range(v.shape[1]):
                sum.append(np.dot(G_i.T, v[:,vi:vi+1]).flatten())
            local_sums.append(sum)
        local_sums = np.asarray(local_sums)
        local_sums = np.sum(local_sums, axis=0).flatten()

        aps = []
        for v, g in zip(Vk_list, G_list_n):
            sum = np.zeros((v.shape[0],1))
            for vi in range(v.shape[1]):
                it = local_sums[vi] * v[:,vi:vi+1].T
                it = np.reshape(it, sum.shape)
                sum = sum + it
            ap = g -sum
            aps.append(ap)

        # concatenate the locally corrected vectors
        G_i = np.concatenate(aps, axis=0)
        G_i = G_i/np.linalg.norm(G_i)
        G_list = []
        start = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1],:])
            start = start+d.shape[1]
        ra, sum = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
    print(iterations)
    G_i = np.concatenate(G_list)
    return G_i
def residuals(V, a=None, sums=None):
    '''
    '''
    if a is None:
        a = sh.generate_random_gaussian(1, V.shape[0])
    sum = np.zeros(V.shape[0])
    if sums is not None:
        for v in range(V.shape[1]):
            sum = sum + s[v] * V[:,v].T
    else:
        for v in range(V.shape[1]):
            sum = sum + np.dot(a, V[:, v:v + 1]) * V[:, v].T
    ap = a - sum
    return ap
def all_eigenvalues_residuals_local(data_list, k=10):
    ug = simulate_guo_1(data_list, maxit=1000)
    u_all = ug
    for i in range(1,k):
        ug2 = simulate_guo_k_residuals_local(data_list, maxit=1000, V_k=u_all)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all


def simulate_guo_k_residuals(local_data, maxit, V_k):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = [] # this are the parital eigevenctors
    iterations = 0
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    k = V_k.shape[1]
    G_i = residuals(V_k).T
    G_i = G_i / np.linalg.norm(G_i)
    Vk_list = []
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        Vk_list.append(V_k[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)
    momentum = 5
    while not ra and iterations<maxit:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], 1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        local_sums = []
        for d, g, v in zip(local_data, G_list, Vk_list):
            G_i = momentum * np.dot(d.T, H_i) + g
            G_list_n.append(G_i)

        # residuals can be computed in a federated manner.
        G_i = np.concatenate(G_list_n, axis=0)
        G_i = residuals(V_k, G_i.T).T
        G_i = G_i/np.linalg.norm(G_i)

        G_list = []
        start = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1],  :])
            start = start+d.shape[1]
        ra, sum = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i
        if momentum >1:
            momentum = momentum - 0.5
    print(iterations)
    G_i = np.concatenate(G_list)
    return G_i

def simulate_guo_k_residuals_local_delayed_ortho(local_data, maxit, V_k):
    '''
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.
    Here we test, whether it is possible to lag 1 round behind with the orthonormalisation
    to avoid an extra communciation round when computing the residuals.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''
    G_list = [] # this are the parital eigevenctors
    iterations = 0
    ra = False
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    k = V_k.shape[1]
    G_i = residuals(V_k).T
    Vk_list = []
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        Vk_list.append(V_k[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not ra and iterations<maxit:
        iterations = iterations + 1

        H_i = np.zeros((local_data[0].shape[0], 1))
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local

        G_i = np.concatenate(G_list, axis=0)
        G_i = residuals(V_k, G_i.T)
        G_i = G_i/np.linalg.norm(G_i)
        G_i = G_i.T

        G_list = []
        start = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]
        G_list_n = []
        for d, g, v in zip(local_data, G_list, Vk_list):
            G_i = np.dot(d.T, H_i) + g
            G_list_n.append(G_i)
        G_list = G_list_n
        ra, sum = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i

    print(iterations)
    G_i = np.concatenate(G_list)
    return G_i

def all_eigenvalues_residuals(data_list):
    ug = simulate_guo_1(data_list, maxit=1000)
    u_all = ug
    for i in range(1,4):
        ug2 = simulate_guo_k_residuals(data_list, maxit=1000, V_k=u_all)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all
def all_eigenvalues_delayed_ortho(data_list):
    ug = simulate_guo_1(data_list, maxit=1000)
    u_all = ug
    for i in range(1,10):
        ug2 = simulate_guo_k_residuals_local_delayed_ortho(data_list, maxit=1000, V_k=u_all)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all

def assure_consecutive(arr):
    '''
    >>> assure_consecutive([0,1,2,3,4,5,6])
    6
    >>> assure_consecutive([1,2,3,4,5,6,8])
    5

    Find the last index such that all numbers
            at previous indices are consecutive natural numbers
    Args:
        arr:

    Returns: the last index

    '''
    if len(arr)==0:
        return -1
    i = 0
    while i < (len(arr) - 1) and arr[i] + 1 == arr[i + 1]:
        i = i + 1
    return i
def simulate_guo_local_rounds(local_data, k, maxit, scipy):
    '''
    , filename, scipy, choices, precomputed_pca=None
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
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    # send parts to local sites
    for d in local_data:
        G_list.append(G_i[start:start + d.shape[1], :])
        start = start + d.shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    while not ra and iterations<maxit:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], k))

        for lr in range(0):
            G_list_n = []
            for d, g in zip(local_data, G_list):
                H_local = np.dot(d, g)
                G_i = np.dot(d.T, H_local) +  1* g
                G_i = G_i/np.linalg.norm(G_i)
                G_list_n.append(G_i)
            G_list = G_list_n
        print(iterations)
        for d, g in zip(local_data, G_list):
            H_local  = np.dot(d, g)
            H_i = H_i + H_local
        G_list_n = []
        for d, g in zip(local_data, G_list):
            G_i = 2* np.dot(d.T, H_i) +  g
            #G_i, Q = la.qr(G_i, mode='economic')
            G_list_n.append(G_i)
        start = 0
        G_i = np.concatenate(G_list_n, axis=0)
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        eigenvals = np.sqrt(eigenvals)
        G_i, R = la.qr(G_i, mode='economic')
        print(co.compute_angles(scipy, G_i))
        G_list = []
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]

        ra, conv = gv.convergence_checker(H_i, H_i_prev)
        H_i_prev = H_i

    print(iterations)
    G_i = np.concatenate(G_list)
    return G_i, eigenvals
import sys

def better_hybrid_scheme(local_data, k, maxit, perc_conv=0.75):
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
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')
    # send parts to local sites
    transmission_logger = []
    di = 0
    for d in local_data:
        transmission_logger.append(
            ["G_i=SC", iterations, di, sys.getsizeof(G_i[start:start + d.shape[1], :].tobytes())])
        G_list.append(G_i[start:start + d.shape[1], :])
        start = start + d.shape[1]
        di = di + 1
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    converged_eigenvals = []
    G_conv = None
    convii = 0
    delta_prev = None
    annormaly_counter = [0]*k
    ann = [0] *k
    nr_conv = np.floor(k*perc_conv)
    eigenvalues_conv = []
    while not ra and iterations < maxit and ann[convii]<10:
        iterations = iterations + 1
        print(iterations)
        print(convii)
        H_i = np.zeros((local_data[0].shape[0], k-convii))
        di = 0
        for d, g in zip(local_data, G_list):
            H_local = np.dot(d, g)
            transmission_logger.append(["H_local=CS", iterations, di, sys.getsizeof(H_local.tobytes())])
            H_i = H_i + H_local
            transmission_logger.append(["H_global=SC", iterations, di, sys.getsizeof(H_i.tobytes())])
            di = di + 1
        G_list_n = []
        ra, conv, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)

        if delta_prev is not None and len(delta) == len(delta_prev):
            for d in range(len(delta)):
                if delta_prev[d]<delta[d]:
                    annormaly_counter[d+convii] = annormaly_counter[d+convii]+1

        delta_prev = delta
        di = 0
        for d, g in zip(local_data, G_list):
            G_i = np.dot(d.T, H_i) + g
            transmission_logger.append(["Gi_local=CS", iterations, di, sys.getsizeof(G_i.tobytes())])
            G_list_n.append(G_i)
            di = di + 1
        start = 0


        G_i = np.concatenate(G_list_n, axis=0)
        eigenvals = eigenvalues_conv.copy()
        diff = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.sqrt(np.linalg.norm(G_i[:, col])))
        for e in range(1,len(eigenvals)):
                diff.append(np.log(np.abs(eigenvals[e - 1] - eigenvals[e])))
        mm = np.mean(diff)
        ssd = np.sqrt(np.var(diff))

        for d in range(len(diff)):
            if diff[d] < mm - ssd:
                ann[d] = ann[d] + 1

        print(ann)

        if G_conv is not None and G_conv.shape[1]>0:
            G_i = np.concatenate([G_conv, G_i], axis = 1)



        G_i, R = la.qr(G_i, mode='economic')

        if len(converged_eigenvals)>0:
            convii = convii+len(converged_eigenvals)
            H_i_prev = H_i[:, len(converged_eigenvals):]
            ann = k*[0]
            eigenvalues_conv = eigenvals[0:convii]

        else:
            H_i_prev = H_i

        G_conv = G_i[:, 0:convii]
        G_i = G_i[:, convii:]

        print(convii)
        print(G_i.shape)
        G_list = []
        di = 0
        for d in local_data:
            G_list.append(G_i[start:start + d.shape[1], :])
            transmission_logger.append(
                ["G_i=SC", iterations, di, sys.getsizeof(G_i[start:start + d.shape[1], :].tobytes())])
            start = start + d.shape[1]
            di = di + 1

        #log_current_accuracy(scipy=scipy, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations,
        #                     filename=filename,
        #                     choices=choices, precomputed_pca=precomputed_pca, transmission_logger=transmission_logger)
        transmission_logger = []
    G_i = np.concatenate(G_list)
    return G_conv, converged_eigenvals

if __name__ == '__main__':
    # parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    # parser.add_argument('-f', metavar='infile', type=str, help='filename of data file; default tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='filename of data file; default tab separated')
    # parser.add_argument('-g', metavar='grm', type=str, default=None)
    # parser.add_argument('-k', metavar='dim', default=10, type=int, help='Number of PCs to calculate')
    # parser.add_argument('-s', metavar='sites', default=10, type=int, help='Number of sites simulated')
    # parser.add_argument('-i', metavar='iteration', default=2000, type=int, help='Maximum number of iterations')
    # parser.add_argument('-p', metavar='outpath', type=str, help='Output directory for result files')
    # args = parser.parse_args()
    # from scipy.sparse import coo_matrix
    # import scaled SNP file
    #data = easy.easy_import(args.f, header=None, rownames=None, center=False, scale_var=False,sep='\t')
    data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    data = coo_matrix.asfptype(data)
    # args.k = 10
   # g = gv.standalone(data, k=2)
    #
    u, s, v = lsa.svds(data.T,k=10)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    data_list, choices = sh.partition_data_vertically(data,2,equal=False, perc=[10,30,60])


    #ug, ev = simulate_guo_local_rounds(data_list, 10, 2000, scipy = u)
    #print(co.compute_angles(u, ug))
    # # produce k-1 eigenvectors
    # start = time.monotonic()
    # u_all = all_eigenvalues_residuals(data_list)
    # print(co.compute_angles(u, u_all))
    # end = time.monotonic()
    # print(end - start)
    #
    #
    # These is the same as above, with the difference
    # that the algorithm to compute the residuals
    # is closer to the implementation required
    # for an actual federated implementation.
    # start = time.monotonic()
    # u_all_loc,con = all_eigenvalues_residuals_local(data_list, 10)
    # print(co.compute_angles(u, u_all_loc))
    # end = time.monotonic()
    # print(end - start)

    start = time.monotonic()
    u_all_loc,con = better_hybrid_scheme(data_list, 10, 2000, perc_conv = 0.8)
    print(co.compute_angles(u, u_all_loc))
    end = time.monotonic()
    print(end - start)
    #
    # # produce k-1 eigenvectors
    # start = time.monotonic()
    # u_all_del = all_eigenvalues_delayed_ortho(data_list)
    # print(co.compute_angles(u, u_all_del))
    # end = time.monotonic()
    # print(end - start)



    # import compress_pickle as cp
    # cpickle = cp.dumps(data, 'gzip')
    # import sys
    # sys.getsizeof(data)/100000000
    # sys.getsizeof(cpickle)/100000000

    from matplotlib import pyplot as plt
    import numpy as np
    import math  # needed for definition of pi

    uu1, ss1, vv1 = lsa.svds(data, k=10)
    uu2, ss2, vv2 = lsa.svds(data_list[0], k=10)
    uu, ss, vv = lsa.svds(data_list[1], k=10)
    np.sqrt(ss)

    print('all')
    ev = np.flip(ss1)
    for i in range(len(ev) - 1):
        print(ev[i]- ev[i + 1])

    print('1')
    ev = np.flip(ss)
    for i in range(len(ev) - 1):
        print(ev[i] -ev[i + 1])

    print('2')
    ev = np.flip(ss2)
    for i in range(len(ev) - 1):
        print(ev[i] - ev[i + 1])
    #
    # x = np.arange(0, 10)
    # plt.plot(x, np.flip(ss))
    # plt.plot(x, np.flip(ss2))
    # plt.plot(x, np.flip(ss1))
    # plt.xlabel("rank")
    # plt.ylabel("eigenvalue")
    # plt.show()
