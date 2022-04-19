import sys
sys.path.append('/home/anne/Documents/featurecloud/pca/federated_dp_pca')
sys.path.append('/home/anne/Documents/featurecloud/pca/scanorama/')
#sys.path.append('/home/anne/Documents/featurecloud/pca/scanorama')
import python.PCA.horizontal as h
import python.PCA.shared_functions as sh
import os.path as op
import os
import python.PCA.horizontal.horizontal_pca_benchmark as hb
import seaborn
import python.PCA.horizontal.horizontal_pca_power_iteration as hori
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
import pandas as pd
from time import time

# ad = 'AD'
# path = '/home/anne/Documents/featurecloud/pca/horizontal-pca/PSO_and_AD_zenodo_upload_revised/PSO_and_AD_zenodo_upload_1sep2020/RNA-Seq/'

ad = ''
path = '/home/anne/Documents/featurecloud/pca/horizontal-pca/PSO_and_AD_zenodo_upload_revised/PSO_and_AD_zenodo_upload_1sep2020/RNA-Seq/psoriasis'

file_list = os.walk(path)
files = [f for f in file_list]

dataframes = []
labels = []
meta_list = []

for d in files:
    for f in d[2]:
        if 'normalized' in f:
            matrix_file = os.path.join(d[0], f)
            da = pd.read_csv(matrix_file, sep='\t', index_col=0)
            try:
                da = da.drop('Unnamed: 17', axis=1)
            except:
                pass
            dataframes.append(da)
            name = op.basename(matrix_file).replace('normalized_counts_matrix_', '').replace('.txt', '')
            labels.append([name]*da.shape[1])
            meta_file = os.path.join(d[0], name+'.txt')
            metadata = pd.read_csv(meta_file, sep='\t', index_col=0)
            metadata = metadata[['GSM', 'diagnosis', 'lesional', 'treatment', 'dose', 'time_point',
                    'time_point_unit', 'tissue', 'anatomical_site', 'platform',
                    'LibraryLayout', 'gender', 'age', 'ethnicity', 'patient_id']]
            meta_list.append(metadata)

labels = np.concatenate(labels)
metadata = np.concatenate(meta_list)
metadata = pd.DataFrame(metadata)
metadata.columns = ['GSM', 'diagnosis', 'lesional', 'treatment', 'dose', 'time_point',
                    'time_point_unit', 'tissue', 'anatomical_site', 'platform',
                    'LibraryLayout', 'gender', 'age', 'ethnicity', 'patient_id']
metadata.to_csv('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/psoriasis/'+ad+'metadata.tsv')



mydata = dataframes[0]
for d in dataframes[1:]:
    mydata = mydata.join(d)

patients = mydata.columns

mydata = mydata.values
mydata = np.log2(mydata+1)
np.isnan(mydata)
means = np.nanmean(mydata, axis=1)
var = np.nanvar(mydata, axis=1)
#highly_var = np.argsort(var)[0:2500]
#mydata = mydata[highly_var, :]
for row in range(mydata.shape[0]):
    mydata[row, np.where(np.isnan(mydata[row,:]))[0]] = means[row]
    mydata[row] =  mydata[row] - means[row]
std = np.nanstd(mydata, axis=1)
for row in range(mydata.shape[0]):
    mydata[row] =  mydata[row] /std[row]

p,c,a = lsa.svds(mydata)
p = np.flip(p, axis=1)
proj = mydata.T @ p

#sns.scatterplot(proj[:,0], proj[:,1], hue=labels, markers=labels)
#plt.show()

proj  = np.concatenate([proj, np.atleast_2d(labels).T, np.atleast_2d(patients).T], axis=1)
proj = pd.DataFrame(proj)
proj.columns = ['PC1', 'PC2','PC3','PC4','PC5','PC6' ,'label', 'patients']
pd.DataFrame(proj).to_csv('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/psoriasis/'+ad+'centralized.tsv')


mydat = []
myproj = []
ind1 = 0
for d in dataframes:
    mydat.append(mydata[:,ind1:(d.shape[1]+ind1)].T)
    myproj.append(proj.iloc[ind1:(d.shape[1]+ind1),:])
    ind1 = d.shape[1]+ind1


myproj = []
ind1 = 0
for d in dataframes:
    sub  = proj.iloc[ind1:(d.shape[1]+ind1),0:6]
    sub = sub.astype(float)
    cov = np.cov(sub.T)
    mean = np.nanmean(sub, axis=0)
    resamples = np.random.multivariate_normal(mean, cov, sub.shape[0])
    myproj.append(resamples)
    ind1 = d.shape[1]+ind1

resamples = np.concatenate(myproj)
resamples  = np.concatenate([resamples, np.atleast_2d(labels).T, np.atleast_2d(patients).T], axis=1)
resamples = pd.DataFrame(resamples)
resamples.columns = ['PC1', 'PC2','PC3','PC4','PC5','PC6' ,'label', 'patients']
pd.DataFrame(resamples).to_csv('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/psoriasis/'+ad+'resampled.tsv')


sns.scatterplot(resamples.iloc[:,0], resamples.iloc[:,1], hue=labels, markers=labels)
plt.show()



u, s = b.simulate_federated_horizontal_pca(mydat)

proj2 = mydata.T @ u
sns.scatterplot(proj2[:,0], proj2[:,1], hue=labels, markers=labels)
plt.show()


proj2  = np.concatenate([proj2[:,0:6], np.atleast_2d(labels).T, np.atleast_2d(patients).T], axis=1)
proj2 = pd.DataFrame(proj2)
proj2.columns =  ['PC1', 'PC2','PC3','PC4','PC5','PC6' ,'label', 'patients']
proj2.to_csv('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/psoriasis/'+ad+'approximate.tsv')


u2, s2, it = hori.simulate_distributed_horizontal(mydat)
proj3 = mydata.T @ u2
sns.scatterplot(proj3[:,0], proj3[:,1], hue=labels, markers=labels)
plt.show()

proj_list = []
for d in mydat:
    u,s,v = lsa.svds(d, k=5)
    proj = d @ np.flip(v.T, axis=1)
    proj_list.append(proj)


projections = np.concatenate(proj_list, axis=0)[:, 0:6]
projections = np.concatenate([projections,np.atleast_2d(labels).T, np.atleast_2d(patients).T], axis=1)
projections = pd.DataFrame(projections)
projections.columns = np.array(['PC1', 'PC2', 'PC3','PC4','PC5','label', 'patient'])
g = sns.FacetGrid(projections, col='label', col_wrap=4)
g.map(sns.scatterplot, 'PC1', 'PC2')
plt.show()

sns.scatterplot(data = projections, x='PC1', y='PC2', hue='label')
plt.show()

pd.DataFrame(projections).to_csv('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/psoriasis/'+ad+'centralized_indiv.tsv')


