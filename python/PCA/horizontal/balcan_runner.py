import balcan as pca
def simulate(datasets):
    partial = []
    for d in datasets:
        partial.append(pca.perform_SVD(d))
    dpca = pca.aggregate_partial_SVDs(partial, 6)
    return dpca