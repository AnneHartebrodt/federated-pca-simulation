import numpy as np
import scipy as sc


def local1(data, G_i):
    return np.dot(data, G_i)

def local2(data, H_pooled, G_previous):
    return G_previous + np.dot(data.T, H_pooled)

def pooling(H_all):
    H_pooled = H_all[0]
    for h in H_all[1:len(H_all)]:
        H_pooled = H_pooled + h
    H_pooled = H_pooled / len(H_all)
    H_pooled = H_pooled / np.linalg.norm(H_pooled)
    return H_pooled

def generate_random_gaussian(n, m):
    draws = n * m
    noise = sc.random.normal(0, 1, draws)
    print('Generated random intial matrix: finished sampling')
    # make a matrix out of the noise
    noise.shape = (n, m)
    # transform matrix into s
    return noise

def standalone(data):
    G_i = generate_random_gaussian(data.shape[1], 1) # phi1
    converged = False
    previous = G_i
    while not converged:
        H_i = np.dot(data, G_i) # YiPhii , gamma if standalone
        G_i =  np.dot(data.T, H_i)  + previous
        G_i = G_i/np.linalg.norm(G_i)
        converged = convergence_checker(G_i, previous)
        previous = G_i
    return G_i

def convergence_checker(current, previous, epsilon=0.0001):
    ra = np.dot(current.T, current)/sc.linalg.norm(current)
    rap = np.dot(previous.T, previous)/sc.linalg.norm(previous)
    if np.abs(ra-rap)<epsilon:
        return True
    else:
        return False