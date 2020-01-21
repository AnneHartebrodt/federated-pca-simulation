import numpy as np
import scipy.linalg as la
import math as math


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
    angle = math.degrees(theta)
    # angle can be Nan
    # print(angle)
    if math.isnan(angle):
        return angle
    # return the canonical angle
    if angle > 90:
        return np.abs(angle - 180)
    else:
        return angle
