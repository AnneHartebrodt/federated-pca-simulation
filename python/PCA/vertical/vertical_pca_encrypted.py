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

import python.PCA.shared_functions as sh
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import python.PCA.comparison as co
import python.PCA.convenience as cv
import time as time
import python.import_export.mnist_import as imnist
from scipy.sparse import coo_matrix
import numpy as np
from Pyfhel import Pyfhel
import tempfile
from pathlib import Path
import python.PCA.vertical.federated_qr as qr
import os
import os.path as op
import pandas as pd

# Using a temporary dir as a "secure channel"
# This can be changed into real communication using other python libraries.
secure_channel = tempfile.TemporaryDirectory()
sec_con = Path(secure_channel.name)
pk_file = sec_con / "mypk.pk"
contx_file = sec_con / "mycontx.con"

##### CLIENT
# HE Object Creation, including the public and private keys
HE = Pyfhel()
HE.contextGen(p=65537, m=2 ** 12)
HE.keyGen()  # Generates both a public and a private key

# Saving only the public key and the context
HE.savepublicKey(pk_file)
HE.saveContext(contx_file)

def log_current_accuracy(u, G_i, eigenvals, conv, current_iteration, filename, choices, precomputed_pca=None,
                         current_ev=None, gi_delta_obj=None, v=None, H_i=None):
    """
    Log the current iterations angle to the canonical
    Args:
        u: column vector matrix with canonical eigenvectors
        G_i: column vector based matrix with current eigenvector estimation
        current_iteration: iteration index
        filename: output filename prefix > out will be saved to out.angles, and out.cor

    Returns: None

    """
    if current_ev is not None:
        info_string = str(current_ev) + '\t' + str(current_iteration)
    else:
        info_string = str(current_iteration)

    if not u is None:
        with open(filename + '.angles.u', 'a+') as handle:
            angles = co.compute_angles(u[choices, :], G_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

        with open(filename + '.cor', 'a+') as handle:
            correlations = co.compute_correlations(u[choices, :], G_i)
            if correlations is not None and len(correlations) > 0:
                info = cv.collapse_array_to_string(correlations, info_string)
                handle.write(info)

    if not v is None:
        with open(filename + '.angles.v', 'a+') as handle:
            angles = co.compute_angles(v, H_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

    with open(filename + '.eigenval', 'a+') as handle:
        info = cv.collapse_array_to_string(eigenvals, info_string)
        handle.write(info)

    if conv is not None:
        with open(filename + '.conv', 'a+') as handle:
            conv = cv.collapse_array_to_string(conv, info_string)
            handle.write(conv)

    if precomputed_pca is not None:
        with open(filename + '.angles_precomp', 'a+') as handle:
            angles = co.compute_angles(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(angles, info_string)
            handle.write(info)

        with open(filename + '.cor_precomp', 'a+') as handle:
            correlations = co.compute_correlations(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(correlations, info_string)
            handle.write(info)
    if gi_delta_obj is not None:
        with open(filename + '.eigenvector_convergence', 'a+') as handle:
            conv = cv.collapse_array_to_string(gi_delta_obj, info_string)
            handle.write(conv)

def log_choices(logfile, filename, choices):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(cv.collapse_array_to_string(choices, filename))


def log_transmission(logfile, log_entry_string, iterations, counter, element, eigenvector=10):

    """
    Dumps object to json to estimate the size.
    Args:
        logfile:
        log_entry_string:
        iterations:
        counter:
        element:
        eigenvector: Dummy value 10 if not put

    Returns:

    """
    with open(logfile + '.transmission', 'a+') as handle:

        if isinstance(element, np.ndarray):
            size = len(element.flatten())
        elif isinstance(element, float) or isinstance(element, int):
            size = 1
        elif isinstance(element, list):
            # cast to numpy array and flatten
            # in case nested list.
            size = len(np.asanyarray(element).flatten())
        else:
            # in case something goes wrong
            size = -1
        handle.write(log_entry_string + '\t' + str(iterations)+ '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(size) + '\n')


def log_time(logfile, algorithm, time, split, repeat):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(algorithm + '\t' + str(split) + '\t' + str(repeat) + '\t' + str(time) + '\n')


def log_costs(filename, action, duration, split, repeat):
    with open(filename+'.encryption', 'a+') as handle:
            handle.write(action+'\t'+str(split)+'\t'+str(repeat)+'\t'+str(duration)+'\n')



####### MATRIX POWER ITERATION SCHEME #######
def simulate_subspace_iteration_encrypted(local_data, k, maxit, filename=None, u=None, choices=None, precomputed_pca=None, fractev=1.0,
                           federated_qr=False, v=None, gradient=True, epsilon=10e-9, log=True, repeat=1):
    """
    Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of dimensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    """




    G_list = []
    iterations = 0
    convergedH = False
    total_len = 0
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')

    # send parts to local sites
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        start = start + local_data[i].shape[1]

    # Initial guess
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)
    G_i_prev = G_i
    converged_eigenvals = []
    shape_H = H_i_prev.shape


    eigenvals_prev = None
    # Convergence can be reached when eigenvectors have converged, the maximal number of
    # iterations is reached or a predetermined number of eignevectors have converged.
    while not convergedH and iterations < maxit and len(converged_eigenvals) < k * fractev:
        iterations = iterations + 1

        # add up the H matrices
        # first client
        H_local = np.dot(local_data[i], G_list[i])
        start = time.monotonic()
        H_i = [HE.encryptFrac(H_local[x, y]) for x in range(H_local.shape[0]) for y in
                       range(H_local.shape[1])]
        end = time.monotonic()
        if log:
            log_costs(filename, 'matrix_encryption', end - start, split=1, repeat=iterations)

        # the other clients
        for i in range(1,len(local_data)):
            # send local H matrices to server
            H_local = np.dot(local_data[i], G_list[i])
            start = time.monotonic()
            H_local_enc = [HE.encryptFrac(H_local[x, y]) for x in range(H_local.shape[0]) for y in
                           range(H_local.shape[1])]
            end = time.monotonic()
            if log:
                log_costs(filename, 'matrix_encryption', end - start, split=1, repeat=iterations)
            start = time.monotonic()
            H_i = [H_i[x] + H_local_enc[x] for x in range(len(H_local_enc))]
            if log:
                log_costs(filename, 'matrix_addition', time.monotonic() - start, split=1, repeat=iterations)


        H_i_dec = [HE.decryptFrac(x) for x in H_i]
        H_i = np.reshape(H_i_dec, shape_H)
        if log:
            log_costs(filename, 'matrix_decryption', time.monotonic() - start, split=iterations, repeat=repeat)

        # free orthonormalisation in terms of transmission cost
        H_i, R = la.qr(H_i, mode='economic')

        # Eigenvector update
        for i in range(len(G_list)):
            # Use gradient based update of the Eigenvectors
            if gradient:
                G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]
            else:
                # Use power iterations based update of the eigenvalue scheme
                G_list[i] = np.dot(local_data[i].T, H_i)

        # This is just for logging purposes
        G_i = np.concatenate(G_list, axis=0)

        # Eigenvalues are the norms of the eigenvecotrs
        eigenvals = []
        for col in range(G_i.shape[1]):
            eigenvals.append(np.linalg.norm(G_i[:, col]))
        # Save eigenvalues
        eigenvals_prev = eigenvals


        G_i, G_list = qr.simulate_federated_qr(G_list, encrypt=True, filename=filename, repeat=repeat, log=log, split=iterations)

        convergedH, deltaH = sh.eigenvector_convergence_checker(H_i, H_i_prev, tolerance=epsilon)
        H_i_prev = H_i
        G_i_prev = G_i
        #if log:
            #log_current_accuracy(u=u, G_i=G_i, eigenvals=eigenvals, conv=deltaH, current_iteration=iterations,
                            # filename=filename, choices=choices, precomputed_pca=precomputed_pca, v=v, H_i=H_i)

    G_i = np.concatenate(G_list)
    print(iterations)
    return G_i


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
    data = data.T
    # args.k = 10
   # g = gv.standalone(data, k=2)
    #
    u, s, v = lsa.svds(data.T,k=2)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    outdir = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption'
    os.makedirs(outdir, exist_ok=True)
    filename = op.join(outdir, 'mnist')

    angles = []
    data_list, choices = sh.partition_data_vertically(data,2)
    uu = simulate_subspace_iteration_encrypted(data_list, k=10, maxit=200, log=True, filename=filename, repeat=1)
    angles.append(co.compute_angles(uu, u))

    angles = np.concatenate(angles, axis=0)
    pd.DataFrame(angles).to_csv(filename+".angles", sep='\t', header=False, index=False)







