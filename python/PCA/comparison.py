import numpy as np
import scipy.linalg as la
import math as math
import python.PCA.convenience as cv
import os.path as path
import scipy.spatial.distance as d
import pandas as pd



def angle(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree. The angle returned
    can maximally be 90, otherwise the complement will be returned
    Args:
        v1: first vector
        v2: second vector

    Returns: angle in degree or NaN

    """
    dp = np.dot(v1, v2)
    norm_x = la.norm(v1)
    norm_y = la.norm(v2)
    co = np.clip(dp / (norm_x * norm_y), -1, 1)
    theta = np.arccos(co)
    a = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(a):
        return a
    # return the canonical angle
    if a > 90:
        return np.abs(a - 180)
    else:
        return a


def compute_angles(canonical, split, reported_angles=20):
    """
    Compute the angles of the vectors in a matrix with matching vectors
    in a second matrix.

    Args:
        canonical: The first matrix of eigenvectors
        split:  The second matrix of eigenvectors
        reported_angles: Limit for number of angles to compute

    Returns: A vector of angles

    """
    angles = list()
    for i in range(min(reported_angles, min(canonical.shape[1], split.shape[1]))):
        a = angle(canonical[:, i], split[:, i])
        angles.append(a)
    return angles

def compute_correlations(canonical, split, reported_angles=20):
    """
        Compute the correlations of the vectors in a matrix with matching vectors
        in a second matrix.

        Args:
            canonical: The first matrix of eigenvectors
            split:  The second matrix of eigenvectors
            reported_angles: Limit for number of correlations to compute

        Returns: A vector of correlations

        """
    correlations = list()
    for i in range(min(reported_angles, min(canonical.shape[1], split.shape[1]))):
        c = np.corrcoef(canonical[:, i], split[:, i])
        correlations.append(c[0,1])
    return correlations


def compute_save_angles(W0, W1, study_id, filename, outfile, reported_angles=20):
    angles = compute_angles(W0, W1, reported_angles=reported_angles)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))

def calculate_euclidean_distance(V1, V2):
    res = []
    for line in range(min(V1.shape[1], V2.shape[1])):
        res.append(d.euclidean(V1[:, line], V2[:, line]))
    return res

def subspace_reconstruction_error_element_wise(data, eigenvectors):
    '''
    returns the average elementwise subspace reconstruction error
    Args:
        data:
        eigenvectors:

    Returns:

    '''
    res = []
    for i in range(eigenvectors.shape[1]):
        proj = np.dot(data, eigenvectors[:, 0:i])
        rec = np.dot(proj, eigenvectors[:, 0:i].T)
        res.append(np.linalg.norm(data - rec) / (data.shape[1] * data.shape[0]))
    return res

def subspace_reconstruction_error(data, eigenvectors):
    '''
    returns the average elementwise subspace reconstruction error
    Args:
        data:
        eigenvectors:

    Returns:

    '''
    res = []
    for i in range(eigenvectors.shape[1]):
        proj = np.dot(data, eigenvectors[:, 0:i])
        rec = np.dot(proj, eigenvectors[:, 0:i].T)
        res.append(np.linalg.norm(data - rec))
    return res

def mev(u, truth):
    k = min(truth.shape[1], u.shape[1])  # number of eigenvectors in subspace
    print(k)
    m = np.dot(u.T, truth)
    sum = 0
    for i in range(k):
        sum = sum + np.linalg.norm(m[:, i], 2)
    mev = sum / k
    return mev

def angle360(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree.
    Args:
        v1: first vector
        v2: second vector

    Returns: angle in degree or NaN

    """
    dp = np.dot(v1, v2)
    norm_x = la.norm(v1)
    norm_y = la.norm(v2)
    co = np.clip(dp / (norm_x * norm_y), -1, 1)
    theta = np.arccos(co)
    a = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(a):
        return a
    # return the canonical angle
    return a