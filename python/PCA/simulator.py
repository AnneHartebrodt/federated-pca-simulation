
import halko_galinsky_runner as hgpca
import guo_vertical_runner as gpca
import numpy as np
import guo_vertical as gv
import shared_functions as sh
import scipy.linalg as la
import import_export.easy_import as easy
import import_export.import_data as imp
import scipy.sparse.linalg as lsa
import comparison as co
import pandas as pd



if __name__ == '__main__':
    data = pd.read_table('/home/anne/Documents/data/radiation/data_new.tsv', sep='\t', header=None)
    datasets = sh.partition_data_horizontally(data.values, 2)


    # standard standalone pca horizontal
    u,s,v = lsa.svds(data, k=6)
    s = np.sqrt(np.flip(s))
    v = np.flip(v.T, axis=1)
    u =np.flip(u, axis=1)


    # balcan pca
    #dpca = bpca.simulate(datasets)
    #co.compute_angles(v, dpca[0])
    #
    #dpcr = pcpca.simulateImtiaz(datasets, seeds.shape[0])
    #co.compute_angles(v, dpcr[0])

    #qupcr = pcpca.simulateQu(datasets, seeds.shape[0])
    #co.compute_angles(v, qupcr[0])

    #k = kpca.simulateKargupta(datasets)
   # co.compute_angles(v, k[0])


    # vertical simulation
    datasetsv = sh.partition_data_vertically(data, 2)

    # standard vertical PCA
    uv,sv,vv = lsa.svds(data.T, k=6)
    sv = np.sqrt(np.flip(sv))
    vv = np.flip(vv.T, axis=1)
    uv =np.flip(uv, axis=1)

    # simulated federated vertical PCA
    uh, sh, vht = hgpca.simulateHalkoGalinskyAdd(datasetsv, p=6, maxit=50)
    hguo, ll = gpca.simulate_guo(datasetsv, k= 20, maxit=2000)


