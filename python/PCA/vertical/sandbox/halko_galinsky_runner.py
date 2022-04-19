import os.path as path
import numpy as np
import scipy.linalg as la
import halko_galinsky as hgpca
import import_export.easy_import as easy
import import_export.import_data as imp
import shared_functions as sh
import scipy.sparse.linalg as lsa
import comparison as co
import pandas as pd

def run_standalone(data, outfile=None, p=10, header=None, rownames=None, center=False, scale_var=False, scale01=False,
                   scale_unit=False, transpose=False, sep='\t', drop_samples=[], log=False, nr_iter=10, qr=True):
    """
    This function performs a regular principal component analysis and saves the result to files containing
    the projection the
    :param datafile: Unscaled datafile
    :param ddpca:
    :param outfile: path and name of the output file without extension
    :param dims: Number of dimensions to return (#eigenvectors and corresponding eigenvalues)
    :param seed: random seed
    :param nr_samples: #variables to select if not all the data columns are to be used for the pca
    :param header: row number which contains the header/ number of header rows
    :param rownames: column number which contains the rownames/sample ids
    :return: projection, eigenvectors and eigenvalues
    """
    # if data is a string, data needs to be read first, otherwise it is
    # assumed to be scaled and ready to use
    if isinstance(data, str):
        data = easy.easy_import(data, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,
                                scale01=scale01, scale_unit=scale_unit, drop_samples=drop_samples, log=log,
                                outfile=outfile, transpose=transpose)

    UT, ST, VT = hgpca.power_method(data=data, p=p, maxit=nr_iter, qr=qr)
    # return the projected datapoints, the eigenvectors and the eigenvalues
    return UT, ST, VT


def simulateHalkoGalinsky(local, p=10, maxit=10):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """

    # order eigenvectors from largest to smallest, to achive
    # ordered eigenvectors

    locals = []
    for l in local:
        locals.append(hgpca.local_step(data=l, maxit=maxit, p=p))
    # Compute the SVD after all the matrices have been computed
    U = hgpca.pooling_step(locals)
    loc_ev = []
    for l in local:
        loc_ev.append(hgpca.local_step_2(data=l, U=U, p=p - 2))

    return loc_ev


def simulateHalkoGalinskyAdd(data_list, maxit, p):
    H_i_list = []
    G_i_list = []
    # Generate inital noise matrices
    total_len = 0
    for d in data_list:
        total_len = total_len +d.shape[1]
    start =  0
    G_i = sh.generate_random_gaussian(total_len, p)
    G_i, R = la.qr(G_i, mode='economic')
    for d in data_list:
        G_i_list.append(G_i[start:start+d.shape[1], :])
        start = start+d.shape[1]
    for i in range(maxit):
        H_i = np.zeros((data_list[0].shape[0], p))
        for d, G_i in zip(data_list, G_i_list): # local1 + global1
            H_i = H_i + np.dot(d, G_i)
        for d in range(len(data_list)): #local2
            #G_i_list[d] = hgpca.compute_Gi(data_list[d], H_i)
            G_i_list[d] = G_i_list[d] + hgpca.compute_Gi(data_list[d], H_i)
        G_i = np.concatenate(G_i_list, axis=0)
        G_i, R = la.qr(G_i, mode='economic')
        G_list = []
        for d in data_list:
            G_i_list.append(G_i[start:start + d.shape[1], :])
            start = start + d.shape[1]
        H_i_list.append(H_i)
    U = hgpca.pooling_step(H_i_list, p - 1) # global2
    result = []
    for d in data_list:
        #result.append(np.dot(U.T, d)) # local 3
        result.append(hgpca.local_step_2(U, d, p-2))
    #full = np.concatenate(result, axis=1)
    #UH, SH, VHT = lsa.svds(full, k=p-2)
    #return UH, SH, VHT
    eigenvalues, l = reformat_SVD_object(result)
    mean = np.mean(eigenvalues, axis=0)
    scale = np.linalg.norm(np.concatenate(l['0']))
    eq = hgpca.eigenvector_assembler_la(data_list, l['0'], scale*mean[0], scale)
    return result, eq, l


def simulateHalkoGalinskySVD(local, p=8, maxit=10):
    """
    Simulate a distributed subspace iteration on a list of
    covariance matrices
    Args:
        local: list of covariance matrices
        p: number of eigenvectors
        tolerance: Error tolerance for convergence criterion

    Returns: The eigenvectors, eigenvalues and the number of iterations until convergence

    """

    G_is = []
    for i in range(len(local)):
        # initialise matrices randomly
        d = local[i].shape[1]
        G_is.append(sh.generate_random_gaussian(d, p))
    H_is = []

    maxit = max(maxit, 1)
    for i in range(maxit):
        for l, G_i in zip(local, G_is):
            # calculate H matrices and append them to the list
            H_is.append(np.dot(l, G_i))
        # combine Hi matrices to one and execute SVD after every step
        H_i = hgpca.pooling_step(H_is)
        H_is = []
        # delete G_is and replace by new guess
        G_is = []
        for l in local:
            G_is.append(hgpca.compute_Gi(l, H_i))

    loc_ev = []
    for l in local:
        loc_ev.append(hgpca.local_step_2(data=l, U=H_i))


    return loc_ev

def reformat_SVD_object(local):
    l = {} # reformatted svd object
    eigenvalues = []
    for i in range(len(local)):
        # collect all local eigenvalues in a common list
        # [[eigenvalues site 1], [eigenvalues site 2], [eigenvalues site 3], ...]
        eigenvalues.append(local[i][1])
        #iterate over the eigenvectors and reformat it
        # {
        # 1: [[eigenvector 1 site 1], [eigenvector 1 site 2],[eigenvector 1 site 3]],
        # 2: [[eigenvector 2 site 1], [eigenvector 2 site 2],[eigenvector 2 site 3]]
        # }
        for e in range(local[i][2].shape[1]):
            if str(e) in l:
                l[str(e)].append(local[i][2][:, e])
            else:
                l[str(e)] = []
                l[str(e)].append(local[i][2][:, e])

    eigenvalues = np.asarray(eigenvalues)
    return eigenvalues, l

def assemble(l):
    eigenvectors_assembled = []
    eigenvector_ass2 = []
    for e in l.keys():
        c, c2 = hgpca.eigenvector_assembler_2_dumb(l[e])
        eigenvectors_assembled.append(c)
        eigenvector_ass2.append(c2)

    eigenvectors_assembled = np.asarray(eigenvectors_assembled).T
    eigenvector_ass2 = np.asarray(eigenvector_ass2).T
    return eigenvectors_assembled, eigenvector_ass2

if __name__ == '__main__':
    print('halko galinksy runner')