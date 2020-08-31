import shared_functions as s
import numpy as np

def aggregate_partial_SVDs(partial, total_dim):
    length = len(partial)
    mymatrix = []
    for l in range(length):
        mymatrix.append([l, np.zeros(l.shape[0], total_dim)])
    blockwise = np.block(partial)
    return blockwise