import python.PCA.shared_functions as sh
import scipy.sparse.linalg as lsa
import python.PCA.comparison as co
import python.PCA.convenience as cv
import time as time
import python.import_export.mnist_import as imnist
from scipy.sparse import coo_matrix
import os
import numpy as np

###################################################################
# This is just a quick sanity check to see if we can get away with
# just running PCA locally.
# Looking a the angles this is clearly not the case.
###################################################################

def central_pca_vs_site_pca(data, splits, file, outdir, k=10):
    """
    This is just to show that the result of the PCA on each site is not the same
    as centralised PCA. Some files are written with some informative measures
    on the accuracy of the PCA
    Args:
        data: Some data set
        splits: The number of sites to simulate
        file: file prefix
        outdir: Output directory
        k: Number of Principal components to retrieve

    Returns: Nothing.

    """

    # Canonical single site PCA
    u, s, v = lsa.svds(data.T, k=k)
    u = np.flip(u, axis=1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    start = time.monotonic()
    with open(os.path.join(outdir, str(start) + '.' + file + '.angles'), 'w') as handle, \
            open(os.path.join(outdir, str(start) + '.' + file + '.eigenvalues'), 'w') as handle2, \
            open(os.path.join(outdir, str(start) + '.' + file + '.correlation'), 'w') as handle3, \
            open(os.path.join(outdir, str(start) + '.' + file + '.subspace'), 'w') as handle4:

        # split the data
        for s in splits:
            data_list, choice = sh.partition_data_vertically(data, s, randomize=True)
            for i in range(20):
                u_d_list = []
                d_d_list = []
                for d in data_list:
                    u_d, s_d, v_d = lsa.svds(d.T, k=k)
                    u_d = np.flip(u_d, axis=1)
                    s_d = np.flip(s_d)
                    v_d = np.flip(v_d.T, axis=1)
                    u_d_list.append(u_d)
                    d_d_list.append(np.dot(d, u_d))
                u_d = np.concatenate(u_d_list, axis=0)

                logstring = str(s) + '\t' + str(i)
                handle2.write(cv.collapse_array_to_string(np.abs(s - s_d), logstring))
                handle.write(cv.collapse_array_to_string(co.compute_angles(u[choice,:], u_d), logstring))
                handle3.write(cv.collapse_array_to_string(co.compute_correlations(u[choice,:], u_d), logstring))


                subs_fed = co.subspace_reconstruction_error(data, u_d)
                logstring = str(s) + '\t' + str(i) + '\tfederated'
                handle4.write(cv.collapse_array_to_string(subs_fed, logstring))

                subs_all = co.subspace_reconstruction_error(data, u)
                logstring = str(s) + '\t' + str(i) + '\tcentralised'
                handle4.write(cv.collapse_array_to_string(subs_all, logstring))


if __name__ == '__main__':
    data, test_lables = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw','train')
    data = coo_matrix.asfptype(data)
    central_pca_vs_site_pca(data, splits=[2], file = 'mnist', outdir='/home/anne/Documents/featurecloud/pca/vertical-pca/results/local_vs_global/mnist/')
