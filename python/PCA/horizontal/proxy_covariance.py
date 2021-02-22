import numpy as np
import python.PCA.horizontal.balcan as b


def aggregate_partial_SVDs(svd_list, t1=10, weights=None,means =None, totaln=None):
    """
    This function aggregates the local proxy covariances by averaging them

    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:

    """

    # Average covariance matrices
    s = len(svd_list)
    # by defaul we take all dimensions available
    if weights is not None:
        Ac = weights[0] * (np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :]))
        for svd in range(1, len(svd_list)):
            Ac = Ac + weights[svd] * (
                np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :]))
    else:
        Ac = np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :])
        for svd in range(1, len(svd_list)):
            Ac = Ac + np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :])
        Ac = 1 / s * Ac

    if means is not None:
        Ac  = Ac+ cross_site_variance(means, weights, totaln)

    U, S, UT, nd = b.svd_sub(Ac, t1)
    S = np.sqrt(S)
    UT = UT.T
    return UT[:, 0:nd], S[0:nd]

def cross_site_variance(means, weights, totaln):
    sum = np.sum(means, axis=0)
    sum = sum/totaln

    mm = np.reshape(means[0]/(totaln*weights[0]) - sum, (len(means[0]),1))
    matr =  np.dot(mm, mm.T) #*weights[0]
    for m in range(len(means)):
        mm = np.reshape(means[m]/(totaln*weights[m]) - sum, (len(means[m]),1))
        #matr = matr + (weights[m]* np.dot(mm, mm.T))
        matr = matr + (np.dot(mm, mm.T))
    return matr

def aggregate_partial_SVDs_means(svd_list, t1=10, weights=None, means = None):
    """
    This function aggregates the local proxy covariances by averaging them

    Function assumes equally shaped covariances matrices.
    :param svd_list: List of local P matrices
    :return:

    """

    # Average covariance matrices
    s = len(svd_list)
    # by defaul we take all dimensions available
    if weights is not None:
        Ac = weights[0] * (np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :]))
        for svd in range(1, len(svd_list)):
            Ac = Ac + weights[svd] * (
                np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :]))
    else:
        Ac = np.dot(svd_list[0][0:t1, :].T, svd_list[0][0:t1, :])
        for svd in range(1, len(svd_list)):
            Ac = Ac + np.dot(svd_list[svd][0:t1, :].T, svd_list[svd][0:t1, :])
        Ac = 1 / s * Ac
    U, S, UT, nd = b.svd_sub(Ac, t1)
    S = np.sqrt(S)
    UT = UT.T
    return UT[:, 0:nd], S[0:nd]








