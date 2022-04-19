import python.PCA.shared_functions as sh
import time
import numpy as np
import scipy.linalg as la
from itertools import chain
import scipy.optimize as optim
import python.import_export.spreadsheet_import as si
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import datasets
from sklearn.decomposition import PCA

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
    delta_list=[]
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
        converged, deltas = sh.eigenvector_convergence_checker(X, X_prev, tolerance=epsilon)
        X2_list.append(X)
        delta_list.append(deltas)
        X_prev = X
    print('converged: ' + str(count))
    return X, E, count, X_list, X2_list, delta_list

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

def leek_the_sneaky_leaker_with_loop(xl1, xl2, d,  dim, delta_list, tolerance=1e-9):
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
        # only append if vector not converged.
        #if delta_list[i][j]<1-tolerance:
        a.append(xl2[i][:, j])
        b.append(xl1[i][d, j])
        total_vars = total_vars + 1
        i = i+1

    a = np.stack(a, axis=0)
    a = a[0:dim, :]
    b = np.array(b[0:dim])
    p, res, rnk, sx = np.linalg.lstsq(a, b)
    return p

def simulate_attack(nfeatures = 50, nsamples=1000, k=10, maxit=500):
    #number of features
    dd = nfeatures
    # number of eigenvectors
    k = k
    maxit=maxit
    print(f'simualting data with '+str(dd)+ ' variables')
    # mat = np.random.random((nsamples, dd))
    # mat2 = np.random.random((nsamples, dd))
    # l = [mat, mat2]
    # #l = sh.scale_datasets(l)
    data, labels = datasets.make_blobs(n_samples=nsamples, n_features=nfeatures)
    l, lab = sh.partition_data_horizontally(data, 2)

    # modified simulated subspace iteration (horizontal)
    print(f'simualting PCA with '+str(k)+ ' eigenvectors and maximally '+str(maxit) + ' iterations')
    v, e, coutn, xlist, x2list, delta_list = simulate_distributed_horizontal(l, maxit=800, p=k)


    ind = int(np.ceil(dd / k))
    #print(f'restricting variable usage to d/k ' + str(ind) + ' iterations')
    #xlist = xlist[0:ind]
    #x2list = x2list[0:ind]

    start = time.monotonic()
    ll = []
    for d in range(dd):
        x= leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd, delta_list=delta_list)
        ll.append(x)
    c1 = np.stack(ll, axis=0)

    # compute the centralized covariance
    cov = np.dot(np.concatenate(l).T, np.concatenate(l))

    pcov = np.dot(np.dot(v, np.diag(e)), v.T)

    acc = np.corrcoef(c1.flatten(), cov.flatten())[0,1]
    pacc = np.corrcoef(pcov.flatten(), cov.flatten())[0,1]

    aer = np.abs(np.sum(np.abs(c1)-np.abs(cov)))
    paer = np.abs(np.sum(np.abs(pcov)-np.abs(cov)))
    print(f'correlation of centralized and reconstructed covariance '+
          str(acc))
    end = time.monotonic()
    elapsed = end - start
    print(f'elapsed time: '+ str(elapsed))

    return acc, pacc, aer, paer, elapsed, delta_list

def run_sample():
    # import some data to play with
    iris = datasets.load_iris()
    diabetes = datasets.load_diabetes()
    breast_cancer = datasets.load_breast_cancer()

    for data in [diabetes, breast_cancer]:
        data = data.data
        print(data.shape)
        dd = data.shape[1]
        data_list= sh.partition_data_horizontally(data, 2)[0]
        l = sh.scale_datasets(data_list)
        k = 2
        # modified simulated subspace iteration (horizontal)
        # print(f'simualting PCA with ' + str(k) + ' eigenvectors and maximally ' + str(maxit) + ' iterations')
        v, e, coutn, xlist, x2list, delta_list = simulate_distributed_horizontal(l, maxit=800, p=k)

        start = time.monotonic()
        ll = []
        for d in range(dd):
            x = leek_the_sneaky_leaker_with_loop(xlist, x2list, d=d, dim=dd, delta_list=delta_list)
            ll.append(x)
        end = time.monotonic()
        c1 = np.stack(ll, axis=0)

        # compute the centralized covariance
        cov = np.dot(np.concatenate(l).T, np.concatenate(l))
        acc = np.corrcoef(c1.flatten(), cov.flatten())[0, 1]
        print(acc)
        print(end-start)


if __name__ == '__main__':
    acc = []
    times = []
    k=10
    counter = 0

    sample=False
    if sample:
        run_sample()
    else:
        with open('/home/anne/Documents/featurecloud/pca/approximative-vertical/results/leek_multi.tsv', 'w') as handle:
            for k in [2,5]:
                for dd in [25, 200, 500,1000,2500]:
                    for s in [1000, 5000, 10000, 20000]:
                        for m, na in zip([500, int(np.floor(dd/(2*k)))], ['unrestricted']):
                            for i in range(10):
                                a,pa,aer, paer, t, dlist = simulate_attack(dd, s, k, m)
                                handle.write(str(counter)+'\t'+str(k)+'\t'+str(dd)+
                                         '\t'+str(s)+'\t'+str(na)+'\t'+str(m)+
                                             '\t'+'reconstructed'+'\t'+str(a)+'\t'+
                                         str(t)+'\n')
                                handle.write(str(counter)+'\t'+str(k)+'\t'+str(dd)+
                                             '\t' + str(s) + '\t' + str(na) + '\t' + str(m) +
                                             '\t' + 'baseline' + '\t' + str(pa) + '\t' +
                                             str(t) + '\n')
                                handle.write(str(counter) + '\t' + str(k) + '\t' + str(dd) +
                                             '\t' + str(s) + '\t' + str(na) + '\t' + str(m) +
                                             '\t' + 'aer' + '\t' + str(aer) + '\t' +
                                             str(t) + '\n')
                                handle.write(str(counter) + '\t' + str(k) + '\t' + str(dd) +
                                             '\t' + str(s) + '\t' + str(na) + '\t' + str(m) +
                                             '\t' + 'paer' + '\t' + str(paer) + '\t' +
                                             str(t) + '\n')


