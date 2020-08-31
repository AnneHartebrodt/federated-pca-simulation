
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



if __name__ == '__main__':
    #seeds = pd.read_table('/home/anne/Documents/data/radiation/data_new.tsv', sep='\t', header=None)
   # seeds = seeds.values[:, 0:7]
    #datasets = sh.partition_data_horizontally(seeds.values, 2)


    # standard standalone pca


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

    file = '/home/anne/Documents/featurecloud/gwas/data/hapmap/thin.scaled.manual'
    outfile = '/home/anne/Documents/featurecloud/gwas/chr10'
    header = 0
    rownames = None
    center = False
    scale = False
    scale_var = False
    scale01 = False
    scale_unit = False
    p = 23
    transpose = True
    sep = '\t'
    #
    import pandas as pd

    data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                            scale01=scale01, scale_unit=scale_unit,
                            outfile=outfile, transpose=transpose, sep=sep)
    data = imp.CustomDataImporter().scale_center_data(data, center=True)
    data = data[:, 0:(data.shape[1] - 1)]
    data = data.T
    datasets = sh.partition_data_vertically(data, 2)


    u,s,v = lsa.svds(data.T, k=6)
    s = np.sqrt(np.flip(s))
    v = np.flip(v.T, axis=1)
    u =np.flip(u, axis=1)

    hg, eq, l = hgpca.simulateHalkoGalinskyAdd(datasets, p=6, maxit=50)
    #v1 = np.flip(l.T)
    #j = 0
    #co.compute_angles(u, v1)
    #co.angle(u[:, 0], np.concatenate([v1[0:5000, 0], -v1[5000:10000, 0]]))
    co.angle(u[:, 0], np.concatenate([l['0'][0], l['0'][1]]))
    co.angle(u[:, 1], np.concatenate([l['1'][0], -l['1'][1]]))

    co.angle(u[:,2], np.concatenate([l['2'][0], l['2'][1]]))
    co.angle(u[:, 3], np.concatenate([l['3'][0], l['3'][1]]))
    hguo, ll = gpca.simulate_guo(datasets, k= 20, maxit=2000)
    co.compute_angles(u, hguo)

    plinkpca = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/hapmap/thin.eigenvec.1', sep='\t')
    plinkpca = plinkpca.values[:, 2:]
    co.compute_angles(plinkpca, hguo)

    import sklearn.decomposition as d
    manual = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/hapmap/thin.scaled.manual', sep='\t',header=None)
    manual = manual.values
    grm = np.dot(manual.T, manual) / manual.shape[0]
    m = d.PCA()
    pca = m.fit(grm)
    pca = pca.components_.T

    co.compute_angles(pca, plinkpca)
    co.compute_angles(plinkpca, hguo)