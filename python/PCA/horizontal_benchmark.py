import python.PCA.power_iteration_runner as por
import python.PCA.shared_functions as sh
import python.import_export.mnist_import as mi
from scipy.sparse import coo_matrix
import numpy as np
import scipy.sparse.linalg as lsa
import python.PCA.comparison as co

if __name__ == '__main__':
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    #
    # data, sample_ids, variable_names = si.data_import('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t', header=0, rownames=0)
    # data = si.scale_center_data_columnwise(data)
    data = coo_matrix.asfptype(data)

    u, s, v = lsa.svds(data, k=10)
    u = np.flip(u, axis=1)
    s = np.flip(s)
    v = np.flip(v.T, axis=1)

    datasplit, choice = sh.partition_data_horizontally(data=data, splits=3)
    x, e, count = por.simulate_distributed(datasplit, p=10, scipy=v)

    print(co.compute_angles(x, v))