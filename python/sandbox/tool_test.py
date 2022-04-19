#import import_export.easy_import as easy
import scipy.sparse.linalg as lsa
import pandas as pd
import python.import_export.mnist_import as imnist
import python.PCA.comparison  as co
from scipy.sparse import coo_matrix

# data = pd.read_table('/home/anne/Documents/data/seeds/seeds_dataset.txt', sep=',')
data, labels = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw/', 'train')
data = coo_matrix.asfptype(data)
import numpy as np
u,s,v = lsa.svds(data)
u  = np.flip(u,1)
spca = pd.read_table('/home/anne/Documents/featurecloud/dev/pca-tool/fed-pca-client/results_2 (1).csv')
spca= spca.values
co.compute_angles(spca[0:60000,],u )