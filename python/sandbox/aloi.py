
import pandas as pd
from scipy.io import arff
import numpy as np
import scipy.sparse.linalg as lsa
import seaborn as sns
import matplotlib.pyplot as plt

from python.PCA.horizontal.leave1out import leave1out_random_sample
def scale_data(data):
    means = np.mean(data, axis = 0)
    #print(means)
    data = data - means
    vars = np.var(data, axis= 0,ddof=1)
    #print(np.sqrt(vars))
    data = data / np.sqrt(vars)
    return data
data = pd.read_csv('/home/anne/Documents/featurecloud/outlier_removal/aloi-27d-50000-max5-tot1508.csv', sep='\s', header=None)

data, meta= arff.loadarff('/home/anne/Documents/featurecloud/outlier_removal/ALOI/ALOI.arff')
da = pd.DataFrame(data)
ol = da.iloc[:,0]
ol = [o.decode('utf-8') for o in ol]
da = da.iloc[:,1:28]
da = da.values

da = scale_data(da)
u, s,v = lsa.svds(da)
v = np.flip(v.T, axis=1)
proj = np.dot(da, v)

proj = pd.DataFrame(proj)
sns.scatterplot(data=proj, x=0, y=1, hue=ol)
plt.legend()
plt.show()

sns.scatterplot(data=proj, x=0, y=2, hue=ol)
plt.show()

sns.scatterplot(data=proj, x=1, y=2, hue=ol)
plt.show()

ao = pd.read_csv('/home/anne/Documents/featurecloud/outlier_removal/ALOI/angle_outliers.tsv', sep ='\t')
aou = ao.values[:,0]
aou = np.unique(aou)
aou = list(aou)

proj['color'] =ol
a = ['a']*proj.shape[0]
for i in range(len(aou)):
    a[aou[i]] = 'ol'
proj['ao'] = a

sns.set_theme(style="ticks")
x = sns.pairplot(proj, hue = 'ao')
x.add_legend()
plt.show()


data_list = [da]
outdir_leave1sample = '/home/anne/Documents/featurecloud/outlier_removal'
name = 'ALOI'
row_list = [range(data.shape[0])]
k = 10
#leave1out_random_sample(data_list, outdir_leave1sample,dataset_name=name, row_list=row_list, k=k)

pd.DataFrame(ol).to_csv('/home/anne/Documents/featurecloud/outlier_removal/ALOI/outliers.tsv')