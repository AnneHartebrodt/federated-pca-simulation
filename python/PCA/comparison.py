import numpy as np
import scipy.linalg as la
import math as math
import convenience as cv
import os.path as path


def angle(v1, v2):
    """
    Calculates the angle between to n-dimensional vectors
    and returns the result in degree
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

def compute_save_angles(W0, W1, study_id, filename, outfile, reported_angles=20):
    angles = compute_angles(W0, W1, reported_angles=reported_angles)
    with open(path.join(outfile, filename), 'a+') as handle:
        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))
