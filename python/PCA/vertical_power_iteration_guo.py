import python.PCA.vertical_pca_library as gv
import python.PCA.shared_functions as sh
import numpy as np
import sys
import numpy as np
import python.PCA.vertical_pca_library as gv
import python.PCA.shared_functions as sh
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import python.import_export.mnist_import as imnist
import python.PCA.comparison  as co
from scipy.sparse import coo_matrix
import time


def residuals(V, a=None, sums=None):
    '''
    '''
    if a is None:
        a = sh.generate_random_gaussian(1, V.shape[0])
    sum = np.zeros(V.shape[0])
    if sums is not None:
        for v in range(V.shape[1]):
            sum = sum + sums[v] * V[:, v].T
    else:
        for v in range(V.shape[1]):
            sum = sum + np.dot(a, V[:, v:v + 1]) * V[:, v].T
    ap = a - sum
    return ap

def simulate_guo(local_data, maxit, V_k=None,starting_vector=None,filename=None, scipy=None, choices=None,  precomputed_pca=None):
    '''
    Retrieve the first eigenvector
    Args:
        local_data: List of numpy arrays containing the data. The data has to be scaled already.
        k: The number of countermensions to retrieve
        maxit: Maximal number of iterations

    Returns: A column vector array containing the global eigenvectors

    '''

    iterations = 0 # iteration counter
    converged = False 
    total_len = 0 # total number of samples/incounterviduals in data
    for d in local_data:
        total_len = total_len + d.shape[1]
    

    if starting_vector is None:
        # if there is no starting vector we generate an orthogonal
        # vector and start iterating
        if V_k is not None:
            # if it is not the first eigenvector use residuals
            # to orthogonalise
            G_i = residuals(V_k).T
        else:
            # if it is the first eigenvecto just generate
            # randonly and scale to unit norm
            G_i = sh.generate_random_gaussian(total_len, 1)
            G_i = G_i / np.linalg.norm(G_i)  # normalize
    else:
        # the vector is already preiterated and is orthonormal, we
        # just have to assure it is a 2d array.
        G_i = np.reshape(starting_vector, (total_len, 1))

    transmission_logger = []  # logger of the size of data to be transmitted.

    start = 0 # running variable to partition G_i
    G_list = []  # this are the parital eigevenctors
    Vk_list = [] # these are the partial eigenvectors already iterated
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        log_transmission(filename, "G_i=SC", iterations, d, G_list[i])
        if V_k is not None:
            Vk_list.append(V_k[start:start + local_data[d].shape[1], :])
        start = start + local_data[d].shape[1]
    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], 1)

    while not converged and iterations < maxit:
        iterations = iterations + 1
        print(iterations)
        H_i = np.zeros((local_data[0].shape[0], 1)) # dummy initialise the H_i matrix
        counter = 0
        for d, g in zip(local_data, G_list):
            H_local = np.dot(d, g)
            log_transmission(filename, "H_local=CS", iterations, counter, H_local)
            H_i = H_i + H_local
            counter = counter+1
        log_transmission(filename, "H_global=SC", iterations, counter, H_i)
        


        for i in range(len(local_data)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]

        gi_norm = 0
        if V_k is None:
            # compute the norm of the eigenvector and done.
            for i in range(len(G_list)):
                local_norm = np.sum(np.square(G_list[i]))
                gi_norm = gi_norm+np.sum(np.square(G_list[i]))
                log_transmission(filename, "local_norm=CS", iterations, i, local_norm)

        else:
            local_sums = []
            for i in range(len(G_list)):
                sum = []
                for vi in range(Vk_list[i].shape[1]):
                    sum.append(np.dot(G_list[i].T, Vk_list[i][:, vi:vi + 1]).flatten())
                local_sums.append(sum)
            local_sums = np.asarray(local_sums)
            local_sums = np.sum(local_sums, axis=0).flatten() # flatten to remove nesting

            for i in range(len(G_list)):
                ap = G_list[i]
                for vi in range(Vk_list[i].shape[1]):
                    it = local_sums[vi] * Vk_list[i][:, vi:vi + 1].T
                    it = np.reshape(it, ap.shape)
                    ap = ap - it
                G_list[i] = ap

            for i in range(len(G_list)):
                c_n = np.sum(np.square(G_list[i]))
                log_transmission(filename, "local_norm=CS", iterations, i, c_n)
                gi_norm = gi_norm + c_n

        gi_norm = np.sqrt(gi_norm)
        for i in range(len(G_list)):
            G_list[i] = G_list[i]/gi_norm

        log_transmission("global_norm=SC", iterations, 1 , gi_norm)

        converged, sum, conv, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)
        H_i_prev = H_i
        log_current_accuracy(scipy, G_i, [gi_norm], conv=delta, current_iteration=iterations,
            filename=filename,choices=choices, precomputed_pca=precomputed_pca, current_ev=1,)

        transmission_logger = [] # reset transmission logger

    # return a complete eigenvector
    G_i = np.concatenate(G_list)
    return G_i

def compute_k_eigenvectors(data_list, k, maxit, filename=None, scipy=None, choices=None, precomputed_pca=None):
    ug = simulate_guo(data_list, maxit=maxit, filename=filename, scipy=scipy, choices=choices)
    u_all = ug
    for i in range(1, k):
        ug2 = simulate_guo(data_list, maxit=maxit, V_k=u_all, filename=filename, scipy=scipy,choices=choices)
        u_all = np.concatenate([u_all, ug2], axis=1)
    return u_all

def log_transmission(logfile, log_entry_string, iterations, counter, element):
    with open(logfile, 'a') as handle:
        handle.write(log_entry_string +'\t'+ str(iterations)
                     + '\t'+ str(counter) +'\t' + str(sys.getsizeof(element.tobytes())))


if __name__ == '__main__':


    data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    data = coo_matrix.asfptype(data)
    # args.k = 10
    # g = gv.standalone(data, k=2)
    #
    u, s, v = lsa.svds(data.T,k=10)
    u = np.flip(u, axis = 1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    data_list, choices = sh.partition_data_vertically(data,2)
    ev = compute_k_eigenvectors(data_list, 10, 2000)
    print(co.compute_angles(u, ev))
