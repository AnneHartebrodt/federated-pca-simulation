import numpy as np
import scipy as sc
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import scipy.spatial.distance as d
import python.PCA.convenience as cv
import time as time


def noise_variance_pow_it(epsilon, p, L, delta):
    '''

    Args:
        epsilon:
        p: iteration rank
        L: number of iterations
        delta:

    Returns: noise variance to be added during each iteration of the noisy power
                method

    '''
    v = np.sqrt(4 * p * L * np.log(1 / delta)) * (1 / epsilon)
    return v


def L_nr_it(eig_k, eig_q1, d):
    '''

    Args:
        eig_k: eigenvalue of target rank
        eig_q1: eigenvalue of intermediate iteration rank
        d: number of dimensions

    Returns: The required numbers of iteration to reach convergence

    '''
    L = np.ceil((eig_k / (eig_k - eig_q1)) * np.log(d))
    return L


def communication_overhead_M(eig_k, eig_q1, d, p, s):
    M = (eig_k / (eig_k - eig_q1)) * s * p * d * np.log(d)
    return M


def coherence(eigenvectors, n):
    '''
    the matrix coherence is defined as follows:
    A=UEU_t (Eigenvalue decomposition)
    A elem mxm
    coh(V) m * ||U||inf
    Args:
        eigenvectors: Eigenvectors from the singular value decomposition as
                        np array
        n: Number of samples

    Returns: The coherence of the matrix

    '''

    m = np.max(eigenvectors)
    coher = n * m
    return coher


def epsilon_bound(coher, d, L, eig_k, eig_q1, epsilon, delta, p):
    '''

    Args:
        coherence: Upper bound on matrix coherence (must be an estimate)
        d: number of dimensions
        L: number of iterations
        eig_k: eigenvalue of target rank
        eig_q1: eigenvalue of intermediate iteration rank
        epsilon:
        delta:
        p: number of dimensions of the data

    Returns: an upper bound for the overall privacy loss after the required
            number of iterations to reach convergence.

    '''
    v = noise_variance_pow_it(epsilon, p, L, delta)
    print(v)
    bound = v * np.sqrt(coher * np.log(d) * np.log(L)) / (eig_k - eig_q1)
    return bound


def generate_random_gaussian(n, m, sigma):
    draws = n * m
    noise = sc.random.normal(0, sigma, draws)
    print('Generated random intial matrix: finished sampling')
    # make a matrix out of the noise
    noise.shape = (n, m)
    # transform matrix into s
    return noise


def variance_explained(eigenvalues, perc=0.5):
    total_variance = sum(eigenvalues)
    percentages = eigenvalues / total_variance
    p = 0
    sum_perc = 0
    while sum_perc < perc and p < len(eigenvalues):
        sum_perc = sum_perc + percentages[p]
        p = p + 1
    return p


def local_step(Xi, data):
    Yi = np.dot(data, Xi)
    return Yi


def pooling_step(Yis, current, weights=None):
    converged = False

    if weights is not None:
        Yi = Yis[0]*weights[0]
        for i in range(1, len(Yis)):
            Yi = Yi + Yis[i]*weights[i]
    else:
        Yi = Yis[0]
        for i in range(1, len(Yis)):
            Yi = Yi + Yis[i]
    # eigenvalues are the column norms of the unnormalised matrix
    E = la.norm(Yi, axis=0)
    Xi, R = la.qr(Yi, mode='economic')

    if np.abs(np.sum(np.dot(np.transpose(Xi), current)) - Xi.shape[1]) < 0.01:
        converged = True
    return Xi, E, converged


def power_method(data, sigma, p, noise=False):
    noise_norms = []
    # U,T,UT = lsa.svds(data, 1000)
    X_0 = generate_random_gaussian(data.shape[0], p, sigma)
    X_0, R = la.qr(X_0, mode='economic')
    nr_iterations = 0
    converged = False
    current = X_0
    while not converged:
        print(nr_iterations)
        # U1, E, UT1 = lsa.svds(X_0,1000)
        nr_iterations = nr_iterations + 1
        if noise:
            G = generate_random_gaussian(data.shape[0], p, np.max(X_0) * sigma)
            noise_norms.append(la.norm(G.flatten()))
            X_0, R = la.qr(np.dot(data, X_0[:, 0:p]) + G, mode='economic')
        else:
            start = time.monotonic()
            X_0, R = la.qr(np.dot(data, X_0[:, 0:p]), mode='economic')
        converged = convergence_checker(X_0[:, 0:p], current[:, 0:p])
        current = X_0
    eigenvals = eigenvalues(X_0, data)
    ord = np.argsort(eigenvals)
    X_0 = np.flip(X_0[:, ord], axis=1)
    eigenvals = np.flip(np.sort(eigenvals))
    proj_global = sc.dot(data, X_0)
    return proj_global, X_0[:, 0:p], eigenvals[0:p], nr_iterations, noise_norms

def power_iteration(data, tolerance = 0.000001):
    """

    Args:
        data:

    Returns:

    """
    X_0 = sc.random.normal(0, 1, data.shape[0])
    nr_iterations = 0
    converged = False
    current =X_0
    while not converged:
        # U1, E, UT1 = lsa.svds(X_0,1000)
        nr_iterations = nr_iterations + 1

        X_0 = np.dot(data, X_0)
        X_0 = X_0/la.norm(X_0)
        if np.abs(np.abs(np.dot(X_0, current)) - 1) < tolerance:
            converged=True
        current = X_0
    eigenval = eigenvalue(data, X_0)
    proj_global = sc.dot(data, X_0)
    return proj_global, X_0, eigenval, nr_iterations


def convergence_checker(current, previous, tolerance=0.000001, required=None):
    '''

    Args:
        current: The current eigenvector estimate
        previous: The eigenvector estimate from the previous iteration
        tolerance: The error tolerance for eigenvectors to be equal
        required: optional parameter for the number of eigenvectors required to have converged

    Returns: True if the required numbers of eigenvectors have converged to the given precision, False otherwise

    '''
    nr_converged = 0
    col = 0
    converged = False
    if required is None:
        required = current.shape[1]
    while col < current.shape[1] and not converged:
        # check if the scalar product of the current and the previous eigenvectors
        # is 1, which means the vectors are 'parallel'
        if np.abs(np.sum(np.dot(np.transpose(current[:, col]), previous[:, col])) - 1) < tolerance:
            nr_converged = nr_converged + 1
            print(nr_converged)
        else:
            break
        if nr_converged >= required:
            converged = True
        col = col + 1
    return converged


def hotelling_deflation(data, eigenvector, eigenvalue, isColumn=True):
    data = data - eigenvalue * np.outer(eigenvector.T, eigenvector)
    return data


def eigenvalue(A, v):
    Av = A.dot(v)
    return v.dot(Av)


def eigenvalues(eigenvectors, cov):
    eigenvals = []
    for v in range(eigenvectors.shape[1]):
        eigenvals.append(eigenvalue(cov, eigenvectors[:, v]))
    return eigenvals


def e_upper(eig_k, eig_q, r, d):
    e_upper = (eig_q / eig_k) * min(1 / np.log(eig_k / eig_q), 1 / np.log(r * d))
    return e_upper


def theorem2_2(noise, eig_k, eig_q1, eigenvectors, p, q, r):
    '''
    These are the constraints that need to be fulfilled for theorem 2.2
    at every iteration
    Args:
        noise: The matrix of gaussian noise
        eig_k: The k'th eignevector
        eig_q1: The q+1th eigenvector
        eigenvectors: The matrix of the top q eigenvectors
        p: the iteration rank
        q: the intermediate rank
        r: some fixed constant r

    Returns: true, if the conditions are fulfilled, noise norm and

    '''
    bound = e_upper(eig_k, eig_q1, r, d)
    noise_norm = la.norm(noise.flatten())
    UqG_norm = la.norm(np.dot(eigenvectors, noise))
    o = bound * (eig_k - eig_q1)
    o2 = o * ((np.sqrt(p) - np.sqrt(q - 1)) / (r * np.sqrt(d)))
    return (5 * noise_norm <= o and UqG_norm <= o2, noise_norm, UqG_norm)


def assumed_noise(eig_k, eig_q1, e):
    return e * (eig_k - eig_q1)


def determine_parameters(cov, epsilon, delta, n, var_exp=0.7, k=10, q1=20, p=40, r=1):
    r = 1  # r is 'some parameter'
    d = cov.shape[1]
    U, E, UT = lsa.svds(cov, n)
    E, indz = cv.extract_eigenvals(E)
    U = np.flip(np.delete(U, indz, axis=1), axis=1)
    k, idx = variance_explained(E, var_exp)
    q1 = 2 * k
    p = 4 * k

    coher = coherence(U, d)
    L = L_nr_it(E[k], E[q1], n)
    noise_variance = noise_variance_pow_it(epsilon, p, L, delta)

    # The question is: are those two epsilons supposed to be the same
    # I don't think so. They just have the same name.
    # this is the epsilon for the noise
    eu = e_upper(E[k], E[q1 - 1], r, d)
    an = assumed_noise(E[k], E[q1], eu)

    # this is the utility parameter
    bound = epsilon_bound(coher, d, L, E[k], E[q1], epsilon, delta, p)
    params = {'coher': coher, 'L': L, 'sigma': noise_variance, 'e_upper': eu, 'assumed_noise': an,
              'bound': bound, 'k': k, 'p': p, 'q1': q1}
    return params


if __name__ == '__main__':
    print('distributed power iteration library file')
