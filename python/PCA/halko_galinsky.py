import numpy as np
import scipy as sc
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import scipy.spatial.distance as d
import convenience as cv
import time as time
import shared_functions as sh
import comparison as  co



def local_step(data, maxit, p):
    G_i = sh.generate_random_gaussian(data.shape[1], p)
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

def pooling_step(H_all, qr=True, p=10):
    H_all = np.concatenate(H_all, axis=1)
    if qr:
        Q, S = la.qr(H_all,mode='economic')
        return Q
    else:
        U, S, V = lsa.svds(H_all, k= p)
        return U

def local_step_2(U, data, p=10):
    T = np.dot(U.T, data)
    UT, ST, VT = lsa.svds(T, k=p)
    VT = np.flip(VT.T, axis=1)
    ST = np.flip(ST)
    UT = np.flip(UT, axis=1)
    return UT, ST, VT


def compute_Gi(data, H_i):
    G_i = np.dot(data.T, H_i) / data.shape[0]
    G_i, R = la.qr(G_i, mode='economic')
    return G_i

def eigenvector_assembler(eigenvectors):
    E = []
    for e in eigenvectors:
        a = co.angle360(e.T, np.random.random(len(e)))
        if a < 90:
            print(str(a) + 'notflipped')
            E = np.concatenate([E, e])
        else:
            E = np.concatenate([E, e * -1])
            print(str(a) + 'flipped')
    E = E / np.linalg.norm(E)
    return E

def eigenvector_assembler_2_dumb(eigenvectors):
    E = np.concatenate([eigenvectors[0], eigenvectors[1]])
    E2 = np.concatenate([eigenvectors[0], -1 * eigenvectors[1]])
    return E, E2

def eigenvector_assembler_la(localdata, localeigenvector, eigenvalue, scale):
    n = len(localeigenvector)
    eq = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            eq[i,j] = np.dot(localdata[j][i,:], localeigenvector[j]/scale)
            if i == j:
                eq[i,j] = eq[i,j]-eigenvalue*localeigenvector[j][i]/scale
    return eq



def find_min(eq):
    import scipy.optimize as opt
    import numpy as np
    matrix = eq
    signum = np.asarray([1]*eq.shape[1])


    def objective(signum, matrix):
        obj = np.sum(np.abs(np.dot(matrix, signum)))
        return obj

    def constraint1(signum):
        return signum[0]*signum[0]-1
    def constraint2(signum):
        return signum[1]*signum[1]-1


    con = {'type':'eq', 'fun':constraint1}
    con2 = {'type': 'eq', 'fun': constraint2}
    cons = [con, con2]

    bounds = [(-1.0, 1.0)]*eq.shape[0]
    sol = opt.minimize(objective, signum, args=(matrix), method='SLSQP', bounds=bounds, constraints=cons)
    print(sol)
    return sol

# paper says 10 iterations are enough
def power_method(data, p=20, maxit =10, qr = True):
    '''
    Runs the algorithm as implemented according to the paper of Galinsky et al
    Args:
        data:
        p:
        maxit:

    Returns:

    '''
    G_i = sh.generate_random_gaussian(data.shape[1], p)
    G_i, R = la.qr(G_i,mode='economic')
    H_list = []
    for i in range(maxit):
        H_i = np.dot(data, G_i)
        G_i = np.dot(data.T, H_i)/data.shape[0]
        G_i, R = la.qr(G_i,mode='economic')
        H_list.append(H_i)
    U = pooling_step(H_list, qr=qr, p = p-1)
    UT, ST, VT = local_step_2(U, data, p = p-2)
    return UT, ST, VT

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

if __name__ == '__main__':
    print('distributed power iteration library file')
