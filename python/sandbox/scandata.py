#%%

import sys
sys.path.append('/home/anne/Documents/featurecloud/pca/federated_dp_pca')
sys.path.append('/home/anne/Documents/featurecloud/pca/scanorama/bin')
#sys.path.append('/home/anne/Documents/featurecloud/pca/scanorama')
import python.PCA.horizontal as h
import python.PCA.shared_functions as sh
import os.path as op
import os
import python.PCA.horizontal.horizontal_pca_benchmark as hb
import seaborn
import python.PCA.horizontal.horizontal_pca_power_iteration as h
import python.PCA.horizontal.balcan as b
import python.PCA.horizontal.bai as bai
import python.PCA.horizontal.proxy_covariance as proxy
import python.PCA.vertical.simulate_federated_vertically_partionned_pca as vertical
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import os.path as op
import python.PCA.comparison as co
import scanorama as scan

from time import time
import numpy as np

NAMESPACE = 'panorama_integrated'

#%%
from process import load_names
#from config import data_names
from scipy.sparse import csr_matrix

#%%

data_names = []
with open('/home/anne/Documents/featurecloud/pca/scanorama/conf/panorama.txt', 'r') as handle:
    for line in handle:
        line = line.strip()
        data_names.append(op.join('/home/anne/Documents/featurecloud/pca/scanorama/',line ))

print(data_names)
datasets, genes_list, n_cells = load_names(data_names[0:1])
datasets, genes = scan.merge_datasets(datasets, genes_list)

labels = []
names = []
curr_label = 0
for i, a in enumerate(datasets):
    labels += list(np.zeros(a.shape[0]) + curr_label)
    names.append(data_names[i])
    curr_label += 1
labels = np.array(labels, dtype=int)

#%%

d = vstack(datasets)

#%%

u, s, vt = lsa.svds(d, k=10)

#%%

print(vt.T.shape)
d.shape

#%%

proj = np.dot(d, vt.T)

#%%

sns.plot
