import pandas as pd
import scipy.sparse.linalg as sc
import plotnine as nine
from plotnine import ggplot, geom_point, aes, stat_smooth, facet_wrap
import matplotlib.pyplot as plt
import numpy as np
from python.import_export import spreadsheet_import as si

data = pd.read_csv('/home/anne/Documents/featurecloud/data/tcga/data_clean/CPTAC-2/coding_only.tsv', sep='\t', header=0)
data= data.values
data=si.scale_center_data(data)

u, s,v = sc.svds(data.T, k=10)
v= np.flip(v.T, axis=1)

pcs = np.dot(data.T, v)

plt.plot(pcs[:,0], pcs[:, 1], 'o')
plt.show()

plt.plot(pcs[:,2], pcs[:, 1], 'o')
plt.show()

plt.plot(pcs[:,2], pcs[:, 0], 'o')
plt.show()

plt.plot(pcs[:,3], pcs[:, 1], 'o')
plt.show()