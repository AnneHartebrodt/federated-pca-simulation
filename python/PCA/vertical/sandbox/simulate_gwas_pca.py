import numpy as np
import scipy as sc
import halko_galinsky as powerit
import halko_galinsky_runner as runner
import sklearn.decomposition as decomp
import import_export.easy_import as easy


def make_test_intervals(n):
    # hardcoded for now
    # more extreme cases maybe later
    unequal_splits = [[0.1, 0.9], [0.3, 0.7], [0.5, 0.5], [0.2, 0.2, 0.2, 0.2, 0.2], [0.1, 0.1, 0.2, 0.2, 0.4],
                      [0.1, 0.1, 0.1, 0.1, 0.6]]
    # , [0.2375, 0.2375, 0.2375, 0.2375, 0.05]]
    interval_end = list()
    sum = 0
    for i in unequal_splits:
        inter = list()
        for k in i:
            sum = sum + k
            inter.append(np.ceil(n * sum))
        sum = 0
        interval_end.append(inter)
    return interval_end, unequal_splits

def chunk_data( data, interval_end):
    data = data.T
    localdata =[]
    start = 0
    for i in range(len(interval_end)):
        end = int(interval_end[i])
        # slice matrix
        data_sub = data[start:end, :]
        # calculate covariance matrix
        start = int(interval_end[i])
        localdata.append(data_sub.T)
    return localdata

def run_split_study(data, p = 20, maxit = 10, interval_end=None):
    # sklean decomposition as baseline
    pca = decomp.PCA(n_components=p)
    pc = pca.fit(data)

    # standalone version as proof that standalone works
    Vs, Ts, nr_iter = runner.run_standalone(data)

    chunks = chunk_data(data, interval_end)
    chu = runner.simulate_distributed(chunks, p=p, maxit=maxit)
    return chu, pc, Vs, Ts


file = '/home/anne/Documents/featurecloud/gwas/data/1000g10/pca.scaled.in'
outfile = '/home/anne/Documents/featurecloud/gwas/chr10/'
header = 0
rownames = None
center = False
scale = False
scale_var = False
scale01 = False
scale_unit = False
p = 10
transpose = True
#
data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                        scale01=scale01, scale_unit=scale_unit,
                        outfile=outfile, transpose=transpose)
data = importer.s
interval_end, unequal_splits = make_test_intervals(data.shape[1])
print(interval_end[2])
chu, pc, Vs, Ts = run_split_study(data = data,p = 20, maxit=10,interval_end=interval_end[2])

pca = decomp.PCA(n_components=p)
pc = pca.fit_transform(data.T)

u,s,v = sc.sparse.linalg.svds(data.T, k=p)
v = np.flip(v.T, axis=1)
pp = np.dot(data.T, v)

proj = np.concatenate([chu[0][2], chu[1][2]], axis=0)
np.linalg.norm(proj[:,0]- pc[:,0])
import scipy as sc