import shared_functions as shared
import kargupta as kargupta
import balcan as b

def simulateKargupta(datasets):
    partial = []
    for d in datasets:
        u,s,v, nd = shared.svd_sub(d, 10)
        partial.append(shared.projection(d, v.T, 6))
    dpca = kargupta.aggregate_partial_SVDs(partial, 6)
    return dpca
