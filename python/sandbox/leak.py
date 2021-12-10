import python.PCA.shared_functions as sh
import time
import numpy as np
import scipy.linalg as la
from itertools import chain

def simulate_distributed_horizontal(local, p=20, maxit=500, epsilon=1e-9):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of data matrices variable X samples
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """
    d = local[0].shape[1]
    X = sh.generate_random_gaussian(d, p)
    X, R = la.qr(X, mode='economic')
    X_prev = X
    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors
    converged = False
    count = 0
    X_list = []
    X2_list = [X]
    while not converged and count < maxit:
        count = count + 1
        # eigenvalues are the column norms of the unnormalised matrix
        X = np.zeros(X.shape)
        for l in local:
            #  are in colums.
            d1 = np.dot(l, X_prev)
            d1 = np.dot(l.T, d1)
            X = X + d1
        E = la.norm(X, axis=0)
        X_list.append(X)
        X, R = la.qr(X, mode='economic')
        X2_list.append(X)
        converged, deltas = sh.eigenvector_convergence_checker(X, X_prev, tolerance=epsilon)
        X_prev = X
    print('converged: ' + str(count))
    return X, E, count, X_list, X2_list

def leek_the_sneaky_leaker(xl1, xl2, d,  dim):
    a = [[xl2[i][:, j] for j in range(xl2[i].shape[1])] for i in range(len(xl2))]
    a = list(chain.from_iterable(a))

    b = [[xl1[i][d, j] for j in range(xl1[i].shape[1])] for i in range(len(xl1))]
    b = list(chain.from_iterable(b))

    a = np.stack(a, axis=0)
    a = a[0:dim, :]
    b = np.array(b[0:dim])
    x = np.linalg.solve(a, b)
    return x

def leek_the_sneaky_leaker_with_loop(xl1, xl2, d,  dim):
    total_vars = 0
    a = []
    b = []
    i=0
    j=0
    #while total_vars < xl1[0].shape[1]*len(xl1):
    while total_vars < dim:
        i = i % (len(xl2)-1)
        if i==0:
            j = j+1
            j = j % (xl1[i].shape[1]-1)
        a.append(xl2[i][:, j])
        b.append(xl1[i][d, j])
        i = i+1
        total_vars = total_vars + 1
    a = np.stack(a, axis=0)
    a = a[0:dim, :]
    b = np.array(b[0:dim])
    x = np.linalg.lstsq(a, b)
    return x[0]

def simulate_attack():
    #number of features
    dd = 50
    # number of eigenvectors
    k = 10
    maxit=800
    print(f'simualting data with '+str(dd)+ ' variables')
    mat = np.random.random((1000, dd))
    mat2 = np.random.random((1000, dd))
    l = [mat, mat2]

    # modified simulated subspace iteration (horizontal)
    print(f'simualting PCA with '+str(k)+ ' eigenvectors and maximally '+str(maxit) + ' iterations')
    v, e, coutn, xlist, x2list = simulate_distributed_horizontal(l, maxit=800)


    ind = int(np.ceil(dd / k))
    print(f'restricting variable usage to d/k ' + str(ind) + ' iterations')
    xlist = xlist[0:ind]
    x2list = x2list[0:ind]

    start = time.monotonic()
    ll = []
    for d in range(dd):
        x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd)
        ll.append(x)
    c1 = np.stack(ll, axis=0)

    # compute the centralized covariance
    cov = np.dot(np.concatenate(l).T, np.concatenate(l))

    acc = np.nansum(np.abs(c1) - np.abs(cov))
    print(f'absolute difference of centralized and reconstructed covariance '+
          str(acc))
    end = time.monotonic()
    elapsed = end - start
    print(f'elapsed time: '+ str(elapsed))

    return acc, elapsed



if __name__ == '__main__':
    acc = []
    times = []
    for i in range(10):
        a,t = simulate_attack()
        acc.append(a)
        times.append(t)
    print(f'average accuracy: '+str(np.nanmean(acc)))
    print(f'average elapsed time: ' + str(np.nanmean(times)))