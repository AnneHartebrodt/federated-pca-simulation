import scipy as sc
from statsmodels import robust


# https://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
def outlier_removal_mad(data, mulitplier=6, nr_pcs=6):
    flagged_pc = list()
    for column in range(min(nr_pcs, data.shape[1])):
        med = sc.median(data[:, column])
        mad = robust.mad(data[:, column])
        for row in range(data.shape[0]):
            if abs(data[row, column] - med) > (mulitplier * mad):
                flagged_pc.append(row)
    flagged_pc = sc.unique(flagged_pc)
    return flagged_pc
