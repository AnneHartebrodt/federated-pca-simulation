import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *
if __name__ == '__main__':

    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    pd.DataFrame(data).to_csv('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/scaled/mnnist.tsv', header=False, index=False, sep='\t')
