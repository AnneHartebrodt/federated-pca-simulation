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
    angles = list()
    for i in range(min(reported_angles, min(canonical.shape[1], split.shape[1]))):
        a = angle(canonical[:, i], split[:, i])
        angles.append(a)
    return angles

def compute_correlations(canonical, split, reported_angles=20):
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

def compute_angles2(plinkfile, qrfile, scikitfile, outputfile):
    angles = []
    plink = pd.read_csv(plinkfile, sep = '\t', header=0)
    qr = pd.read_csv(qrfile, sep='\t', header=None)
    scikit = pd.read_csv(scikitfile, sep='\t', header=None)

    for i in range(min(min(plink.shape[1], qr.shape[1]), scikit.shape[1])):
        angles.append(['plink vs. qr', str(i), angle(plink.iloc[:, 2+i], qr.iloc[:, i])])
        angles.append(['plink vs. scikit', str(i), angle(plink.iloc[:, 2+i], scikit.iloc[:, i])])
        angles.append(['scikit vs. qr', str(i), angle(scikit.iloc[:, i], qr.iloc[:,i])])
    angles = np.asarray(angles)
    angles = pd.DataFrame(angles)
    angles.to_csv(outputfile, sep='\t', header=None, index=None)