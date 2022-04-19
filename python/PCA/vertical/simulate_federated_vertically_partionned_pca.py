"""
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

"""

##########################################################################################
#
#    If you just want some pseudo code for federated vertically partionned PCA,
#    this here is your script. No unnecessary logging or other stuff.
#
##########################################################################################

import numpy as np
import pandas as pd
import scipy.sparse.linalg as lsa
import scipy.linalg as la
import python.PCA.shared_functions as sh
import python.import_export.mnist_import as mi

import python.PCA.vertical.simulate_federated_qr_orthonormalisation as qr
import python.PCA.comparison as co
import python.import_export.spreadsheet_import as si
import python.PCA.vertical.approximate_vertical_pca_benchmark as approx

def simulate_subspace_iteration(local_data, k, maxit, federated_qr=False,epsilon=10e-9):
    """
    Simulate a federated run of principal component analysis for federated algorithm.

    The input to the function is a list of numpy arrays of the shape [(d, n_1), (d, n_2), ..., (d, n_i)].
    This means that the 'fixed' dimension is on axis 0 and the 'distributed' dimension is on axis 2.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations
        federated_qr: Use federatated QR orthonormalisation

    Returns:
        G_i: the right singular vectors
        eigenvalues: the eigenvalues
        H_i: the left singular vectors

    """

    G_list = []

    iterations = 0
    converged = False

    # generate an intitial  orthogonal noise matrix
    total_len = 0
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0

    # inital orthogonalisation
    # this can be done at one site, since it is random data
    G_i = sh.generate_random_gaussian(total_len, k)
    #G_i, R = la.qr(G_i, mode='economic')

    # send parts to local sites
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        start = start + local_data[i].shape[1]

    # Initial guess (just for comparison)
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)

    # Convergence can be reached when eigenvectors have converged, the maximal number of
    # iterations is reached or a predetermined number of eignevectors have converged.
    while not converged and iterations < maxit:
        print(iterations)
        iterations = iterations + 1
        # add up the H matrices
        # init 0 matrix
        H_i = np.zeros((local_data[0].shape[0], k))
        for i in range(len(local_data)):
            # 'send' local H matrices to server
            H_local = np.dot(local_data[i], G_list[i])
            # add up H matrices at server and send them back to the clients
            H_i = H_i + H_local

        # orthonormalise H
        H_i, R = la.qr(H_i, mode='economic')

        # Compute eigenvalues
        for i in range(len(G_list)):
            # Use power iterations based update of the eigenvalue scheme
            G_list[i] = np.dot(local_data[i].T, H_i)

        # This is just for logging purposes
        G_i = np.concatenate(G_list, axis=0)

        # Eigenvalues are the norms of the eigenvectors
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))

        # if iterations % 10 == 0:
        #     if federated_qr:
        #         G_i, G_list, r, rl = qr.simulate_federated_qr(G_list)
        #     else:
        #         G_i, R = la.qr(G_i, mode='economic')

        # check convergence
        converged, delta = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        # save H to compare at next iteration
        H_i_prev = H_i
    print('converged: '+ str(iterations))
    return G_i, eigenvals, H_i, iterations


####### BENCHMARK RUNNER #######

if __name__ == '__main__':

    path = '~/Documents/featurecloud/test-environment/controller/data/app_test/data/data.tsv'
    #data, test_lables = mi.load_mnist(path, 'train')
    #data = coo_matrix.asfptype(data)

    #data = pd.read_csv(path, sep='\t', header=0, index_col=0).values
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # Transpose data such, that
    data = si.scale_center_data_columnwise(data)
    data = data.T

    u, s, v = lsa.svds(data, k=9)

    v = np.flip(v.T, axis = 1)

    local_data, choice = sh.partition_data_vertically(data=data, splits=3, randomize=True)

    # G, ev, H, it = simulate_subspace_iteration(local_data, k=10, maxit=500, federated_qr=False)
    # co.compute_angles(G, v[choice, :])
    # co.compute_angles(H, np.flip(u, axis=1))
    #
    # G, ev, H, it = simulate_subspace_iteration(local_data, k=10, maxit=500, federated_qr=True)
    # co.compute_angles(G, v[choice, :])
    # co.compute_angles(H, np.flip(u, axis=1))

    gi = approx.run_randomized(local_data, 10,10,500, use_approximate=False, factor_k=2)
    co.compute_angles(gi, v[choice, :])