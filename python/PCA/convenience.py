import os as os
import os.path as path
import time as time

def collapse_array_to_string(a, study_id):
    res = study_id + '\t'
    for e in a:
        res = res + str(e) + '\t'
    res = res + '\n'
    return res


def make_eigenvector_path(inputfile, foldername):
    """
    creates a folder called eigenvectors in the input directory
    if the eigenvector folder already exists a folder named
    'eigenvectors<currentTimestamp> will be created
    Args:
        inputfile: name of the inputfile

    Returns: a pathname according to the stated rules

    """
    print(path.dirname(inputfile))
    if not os.path.exists(path.dirname(inputfile) + '/' + foldername):
        pn = path.dirname(inputfile) + '/' + foldername
        os.makedirs(pn)
    else:
        print('Path exists')
        pn = path.dirname(inputfile) + '/' + foldername + str(time.process_time())
        os.makedirs(pn)
    return pn

def extract_eigenvals(E):
    '''
    Eigendecomposition from scipy.linalg.sparse returns eigenvalues ordered in
    increasing order, followed by eigenvalues which are 0.
    Eigenvalues are returned in decreasing order ommiting th 0s alltogether
    Args:
        E: Eigenvalues from a sparse singular value decomposition.

    Returns: Eigenvalue vector in decreasing order, without 0s.

    '''
    E = E[E != 0]
    return E#, indz