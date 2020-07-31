import numpy as np
import scipy as sc
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import scipy.spatial.distance as d
import convenience as cv
import time as time

def variance_explained(eigenvalues, perc=0.5):
    total_variance = sum(eigenvalues)
    percentages = eigenvalues / total_variance
    p = 0
    sum_perc = 0
    while sum_perc < perc and p < len(eigenvalues):
        sum_perc = sum_perc + percentages[p]
        p = p + 1
    return p



def generate_random_gaussian(n, m):
    draws = n * m
    noise = sc.random.normal(0, 1, draws)
    print('Generated random intial matrix: finished sampling')
    # make a matrix out of the noise
    noise.shape = (n, m)
    # transform matrix into s
    return noise


def local_step(data, maxit, p):
    G_i = generate_random_gaussian(data.shape[1], p)
    G_i, R = la.qr(G_i,mode='economic')
    nr_iterations = 0
    H_list = []
    for i in range(maxit):
        print(nr_iterations)
        nr_iterations = nr_iterations + 1
        H_i = np.dot(data, G_i)
        G_i = np.dot(data.T, H_i) / data.shape[0]
        G_i, R = la.qr(G_i,mode='economic')
        H_list.append(H_i)
    H_all = np.concatenate(H_list, axis=1)
    return H_all


def pooling_step(H_all, p =10):
    H_all = np.concatenate(H_all, axis=1)
    Q, S = la.qr(H_all,mode='economic')
    #U, S, V = lsa.svds(H_all, k=p)

    #U = np.flip(U, axis=1)

    return Q

def local_step_2(U, data, p=5):
    T = np.dot(U.T, data)
    UT, ST, VT = lsa.svds(T, k=p)
    VT = np.flip(VT.T, axis=1)
    ST = np.flip(ST)
    UT = np.flip(UT, axis=1)
    return UT, ST, VT


def local_1(data, G_i):
    H_i = np.dot(data, G_i)
    return H_i

def local_2(H_i, data):
    G_i = np.dot(data.T, H_i) / data.shape[0]
    return G_i




def convergence_checker(current, previous, tolerance=0.000001, required=None):
    '''

    Args:
        current: The current eigenvector estimate
        previous: The eigenvector estimate from the previous iteration
        tolerance: The error tolerance for eigenvectors to be equal
        required: optional parameter for the number of eigenvectors required to have converged

    Returns: True if the required numbers of eigenvectors have converged to the given precision, False otherwise

    '''
    nr_converged = 0
    col = 0
    converged = False
    if required is None:
        required = current.shape[1]
    while col < current.shape[1] and not converged:
        # check if the scalar product of the current and the previous eigenvectors
        # is 1, which means the vectors are 'parallel'
        if np.abs(np.sum(np.dot(np.transpose(current[:, col]), previous[:, col])) - 1) < tolerance:
            nr_converged = nr_converged + 1
            print(nr_converged)
        else:
            break
        if nr_converged >= required:
            converged = True
        col = col + 1
    return converged


# paper says 10 iterations are enough
def power_method(data, p=20, maxit =10):
    G_i = generate_random_gaussian(data.shape[1], p)
    G_i, R = la.qr(G_i,mode='economic')
    nr_iterations = 0
    H_list = []
    for i in range(maxit):
        nr_iterations = nr_iterations + 1
        H_i = np.dot(data, G_i)
        G_i = np.dot(data.T, H_i)/data.shape[0]
        G_i, R = la.qr(G_i,mode='economic')
        H_list.append(H_i)
    H_all = np.concatenate(H_list, axis=1)
    U,S,V = lsa.svds(H_all, k=p-1)
    T = np.dot(U.T, data)
    UT, ST, VT = lsa.svds(T, k=p-2)
    VT = np.flip(VT.T, axis=1)
    ST = np.flip(ST)
    return UT, ST, VT


def power_method_via_qr(data, p=20, maxit =10):
    G_i = generate_random_gaussian(data.shape[1], p)
    G_i, R = la.qr(G_i,mode='economic')
    nr_iterations = 0
    H_list = []
    for i in range(maxit):
        nr_iterations = nr_iterations + 1
        H_i = np.dot(data, G_i)
        G_i = np.dot(data.T, H_i)/data.shape[0]
        G_i, R = la.qr(G_i,mode='economic')
        H_list.append(H_i)
    H_all = np.concatenate(H_list, axis=1)
    Q,S = la.qr(H_all,mode='economic')
    T = np.dot(Q.T, data)
    UT, ST, VT = lsa.svds(T, k=p-2)
    VT = np.flip(VT.T, axis=1)
    ST = np.flip(ST)
    return UT, ST, VT



if __name__ == '__main__':
    print('distributed power iteration library file')
