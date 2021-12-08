import python.PCA.shared_functions as sh
import scipy.sparse.linalg as lsa
import scipy.linalg as la
import scipy as sc
import numpy as np
from python.PCA.logging import *
from python.PCA.vertical.vertical_pca_benchmark import simulate_subspace_iteration

import argparse as ap
import os
import os.path as op
import time
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.sparse.linalg as lsa
from scipy.sparse import coo_matrix
import sys
from itertools import chain
import python.import_export.mnist_import as mi

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
        #log_current_accuracy(v=v, H_i=X, eigenvals=E, conv=None, current_iteration=count, filename=filename,
         #                    choices=choices, precomputed_pca=precomputed_pca, u=None, G_i=None)
        X_prev = X
    print('converged: ' + str(count))
    return X, E, count, X_list, X2_list

def leek_the_sneaky_leaker(xl1, xl2, d,  dim):
    a = [[xl2[i][:, j] for j in range(xl2[i].shape[1])] for i in range(len(xl2))]
    a = list(chain.from_iterable(a))
    print(a)
    b = [[xl1[i][d, j] for j in range(xl1[i].shape[1])] for i in range(len(xl1))]
    b = list(chain.from_iterable(b))
    print(b)
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
    while total_vars < xl1[0].shape[1]*len(xl1):
        i = i % (len(xl2)-1)
        if i==0:
            j = j+1
            j = j % (xl1[i].shape[1]-1)
        a.append(xl2[i][:, j])
        b.append(xl1[i][d, j])
        i = i+1
        total_vars = total_vars + 1
    a = np.stack(a, axis=0)
    print(a.shape)
    a = a[0:dim, :]
    b = np.array(b[0:dim])
    x = np.linalg.lstsq(a, b)
    return x[0]

if __name__ == '__main__':
    mat = np.random.random((90,10))
    mat2 = np.random.random((90,10))
    l = [mat, mat2]

    u, s, vt = la.svd(np.concatenate(l, axis=0))


    v , e, coutn, xlist, x2list = simulate_distributed_horizontal(l, maxit=800)

    dd = 10
    # x = leek_the_sneaky_leaker(xlist, x2list, d=0, dim=dd)
    # x1 = leek_the_sneaky_leaker(xlist, x2list, d=1, dim=dd)
    # x2 = leek_the_sneaky_leaker(xlist, x2list, d=2, dim=dd)
    # c1 = np.stack([x,x1,x2])
    # cov = np.dot(np.concatenate(l).T,np.concatenate(l))
    # xt = x.T

    ll = []
    dd = 10
    for d in range(dd):
        x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd)
        ll.append(x)
    c1 = np.stack(ll, axis=0)
    cc = np.tril(c1)
    cc1 = np.tril(c1, -1).T +np.tril(c1)
    cov = np.dot(np.concatenate(l).T, np.concatenate(l))
    print(np.nansum(np.abs(cc1)-np.abs(cov)))
    print(np.nansum(np.abs(c1) - np.abs(cov)))
    a = c1-cov
    b =cc1-cov
    from sklearn import datasets


    from sklearn.decomposition import PCA

    # import some data to play with
    iris = datasets.load_iris()
    iris = iris.data.T
    v, e, coutn, xlist, x2list = simulate_distributed_horizontal([iris], maxit=800)
    ll = []
    cov = np.dot(iris.T, iris)
    dd = iris.shape[1]
    for d in range(dd):
        x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd)
        ll.append(x)
    c1 = np.stack(ll, axis=0)
    print(np.nansum(c1 - cov))
    #
    x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=2, dim=dd)

    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = np.delete(data, np.where(np.nansum(data, axis=0)==0), axis=1)
    cov = np.dot(data.T, data)
    v, e, coutn, xlist, x2list = simulate_distributed_horizontal([data], maxit=800)
    ll = []
    dd = data.shape[1]
    for d in range(50):
        x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd)
        ll.append(x)
    c1 = np.stack(ll, axis=0)
    np.nansum(c1-cov)

# x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=45, dim=v.shape[0])
    # print(x[0]-cov[45, :])
    # np.nansum(x[0]-cov[45, :])
