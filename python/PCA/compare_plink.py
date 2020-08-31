import os.path as path
import numpy as np
import scipy.linalg as la
import halko_galinsky as hgpca
import import_export.easy_import as easy
import import_export.import_data as imp
import shared_functions as sh
import scipy.sparse.linalg as lsa
import comparison as co
import pandas as pd
import halko_galinsky_runner as hgr


def run_scipy(data, final_dim):
    u, s, v = lsa.svds(data, k=final_dim)
    u = np.flip(u, axis=1)
    v = np.flip(v.T, axis=1)
    s = np.flip(s)
    return u, s, v


def make_result(eigenvectors, v):
    res = []
    for i in range(v.shape[1]):
        res.append(co.angle(eigenvectors[:, i], v[:, i]))
    for i in range(v.shape[1]):
        res.append(co.mev(eigenvectors[:, 0:i + 1], v[:, 0:i + 1]))
    return res

def compare_standalone_to_plink_and_scipy(data, eigenvector_plink, outfile, p, nr_iter, nr_repeats, header=False, rownames=False, center=False, scale_var=False, scale01 = False, scale_unit = False, transpose=False):
    # make a results dataframe containing
    # for both qr and
    # - subspace reconstruction error (1)
    # - mev (final_dim)
    # - angles between eigenvectors (final_dim)
    # - eigenvalues (final_dim)

    eigenvalues_qr = []
    eigenvalues_svd = []

    subspace_rec_err_qr = []
    subspace_rec_err_svd = []

    results_qr = []
    results_svd = []

    results_qr_vs_plink = []
    results_svd_vs_plink = []

    u,s,v = run_scipy(data, final_dim)

    # repeat the standalone PCA to be sure the result is not influenced by random init
    for i in range(nr_repeats):
        UT_qr, ST_qr, VT_qr = hgr.run_standalone(data, outfile=outfile, p=p, header=header, rownames=rownames,
                                             center=center, scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,
                                             transpose=transpose, nr_iter=nr_iter, qr=True)

        UT_svd, ST_svd, VT_svd = hgr.run_standalone(data, outfile=outfile, p=p, header=header, rownames=rownames,
                                                center=center, scale_var=scale_var, scale01=scale01,
                                                scale_unit=scale_unit,
                                                transpose=transpose, nr_iter=nr_iter, qr=False)

        results_qr.append(make_result(VT_qr, v))
        results_svd.append(make_result(VT_svd, v))

        results_qr_vs_plink.append(make_result(VT_qr, eigenvector_plink))
        results_svd_vs_plink.append(make_result(VT_svd,eigenvector_plink))

        eigenvalues_qr.append(np.sqrt(ST_qr))
        eigenvalues_svd.append(np.sqrt(ST_svd))

        subspace_rec_err_qr.append(co.subspace_reconstruction_error(data, VT_qr))
        subspace_rec_err_svd.append(co.subspace_reconstruction_error(data, VT_svd))

    pd.DataFrame(eigenvalues_qr).to_csv(path.join(outfile, 'eigenvalues_qr.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(eigenvalues_svd).to_csv(path.join(outfile, 'eigenvalues_svd.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(subspace_rec_err_qr).to_csv(path.join(outfile, 'sub_rec_err_qr.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(subspace_rec_err_svd).to_csv(path.join(outfile, 'sub_rec_err_svd.tsv'), sep='\t', header=False, index=False)


    pd.DataFrame(np.asarray(results_svd_vs_plink)).to_csv(path.join(outfile, 'svd_vs_plink.tsv'), sep='\t', header=False,index=False)
    pd.DataFrame(np.asarray(results_qr_vs_plink)).to_csv(path.join(outfile, 'qr_vs_plink.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(np.asarray(results_qr)).to_csv(path.join(outfile, 'scipy_vs_qr.tsv'), sep='\t', header=False, index=False)
    pd.DataFrame(np.asarray(results_svd)).to_csv(path.join(outfile, 'scipy_vs_svd.tsv'), sep='\t', header=False, index=False)


if __name__ == '__main__':

    #file = '/home/anne/Documents/data/mnist/flat.csv'
    file = '/home/anne/Documents/data/radiation/data_new.tsv'
    outfile = '/home/anne/Documents/featurecloud/gwas/chr10'
    header = 0
    rownames = None
    center = False
    scale = False
    scale_var = False
    scale01 = False
    scale_unit = False
    p = 12
    transpose = False
    sep = '\t'

    # import and center data, because this is apparently done automaticall in fastPCA
    data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                            scale01=scale01, scale_unit=scale_unit,
                            outfile=outfile, transpose=transpose, sep=sep)
    #data = data[0:200,:]
    data = imp.CustomDataImporter().scale_center_data(data, center=True)


    # eigenvector_plink = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/1000g10/plink2.eigenvec', header=0)
    # eigenvector_plink = np.array(eigenvector_plink)
    # eigenvector_plink = eigenvector_plink[:, 2:eigenvector_plink.shape[1]]

    # eigenvalues_plink = pd.read_table('/home/anne/Documents/featurecloud/gwas/data/1000g10/plink2.eigenval',
    #                                 header=None)
    # eigenvalues_plink = np.asarray(eigenvalues_plink)



    # regular PCA as implemented in scipy
    # as we use the sparse version of scipy we need to flip the order of the eigenvectors
    final_dim = 10

    p = 10
    chunks = [0.5, 1.0]
    result_fed = []
    result_fed_2 = []
    result_fed_vs_plink = []
    eigenvalues_fed = []
    subspace_rec_fed = []
    ut, st, vt =hgr.run_standalone(data)
    u,s,v = run_scipy(data, 3)
    loca = sh.partition_data_vertically(data, 3)
    loc_ev = hgr.simulateHalkoGalinsky(loca, maxit=20, p=p)
    res, eq = hgr.simulateHalkoGalinskyAdd(loca, maxit=30, p=p)
    #Q, R = la.qr(E1)
    #Q2 , R2 = la.qr(E2)
    #loc_ev = hgr.simulateHalkoGalinskySVD(loca, maxit=20, p=p)
    a,b = hgr.reformat_SVD_object(res)

    import itertools as iter


    def brute_force_minimizer(eq):
        com = iter.product([1, -1], repeat=3)
        min = 10000000000
        sol = [0,0]
        for i in com:
            s = np.sum(np.abs(np.dot(eq, np.asarray(i))))
            if s < min:
                print(s)
                print(i)
                sol = np.asarray(i)
                min = s
        return sol
    brute_force_minimizer(eq)
    #hgpca.find_min(eq)
    v1 = np.concatenate([b['0'][0], b['0'][1], b['0'][2]])
    v1 = v1/np.linalg.norm(v1)
    v2 = np.concatenate([b['1'][0], -b['1'][1],b['0'][2]])
    v3 = np.concatenate([b['2'][0], b['2'][1],b['0'][2]])
    co.subspace_reconstruction_error(data, np.asarray([v1,v2, v3]).T)
    co.angle(v1.T, v[:,0])

    np.dot(v1.T, v[:,0])
    np.dot(v2.T, v[:, 0])
    #co.angle(np.concatenate([b['0'][0],b['0'][1], b['0'][2], b['0'][3], b['0'][4]]), v[:, 0])
    # # result_fed.append(make_result(eigenvectors_assembled, v))
    # # result_fed_2.append(make_result(eigenvector_ass2, v))
    # # result_fed_vs_plink.append(make_result(eigenvectors_assembled, eigenvector_plink))
    #
    # subspace_rec_fed.append(co.subspace_reconstruction_error(data, eigenvectors_assembled))
    # Q, R = la.qr(eigenvectors_assembled)
    # Q2, R = la.qr(eigenvector_ass2)
    #
    # pd.DataFrame(eigenvalues_fed).to_csv(path.join(outfile, 'eigenvalues_fed.tsv'), sep='\t', header=None, index=None)
    # pd.DataFrame(subspace_rec_fed).to_csv(path.join(outfile, 'sub_rec_err_fed.tsv'), sep='\t', header=None, index=None)
    # result_fed = np.asarray(result_fed)
    # result_fed_vs_plink = np.asarray(result_fed_vs_plink)
    # pd.DataFrame(result_fed_vs_plink).to_csv(path.join(outfile, 'fed_vs_plink.tsv'), sep='\t', header=None, index=None)
    # pd.DataFrame(result_fed).to_csv(path.join(outfile, 'scipy_vs_fed.tsv'), sep='\t', header=None, index=None)