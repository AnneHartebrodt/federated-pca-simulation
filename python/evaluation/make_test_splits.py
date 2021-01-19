from python.PCA.shared_functions import partition_data_vertically
import pandas as pd
import os
import numpy as np

from python.import_export.mnist_import import load_mnist

def split_mnist():
    dirname = '/home/anne/Documents/featurecloud/dev/pca-tool/fed-pca-client/test_datasets/mnist'
    os.makedirs(dirname, exist_ok=True)
    data, test_lables = load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    data = data.astype('float')
    split_and_write(data, dirname, 'mnist')

def split_compass():
    dirname = '/home/anne/Documents/featurecloud/dev/pca-tool/fed-pca-client/test_datasets/compass'
    os.makedirs(dirname, exist_ok=True)
    data = pd.read_csv('/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv', sep='\t')
    colnames= data.columns
    data = data.values
    split_and_write(data, dirname, 'compass', colnames=colnames)



def split_and_write(data, directory, name, colnames=None):
    data_splits, choices = partition_data_vertically(data, splits=2, equal=True, randomize=False)
    start = 0
    for i in range(len(data_splits)):
        d = pd.DataFrame(data_splits[i])
        names = ['R' + str(x) for x in range(d.shape[0])]
        d.insert(0, 'rownames', names)
        d.to_csv(os.path.join(directory, name) + '_split_' + str(
            i) + '.tsv', header=False, index=False, sep='\t')
        d = d.sample(frac=1)
        d.to_csv(os.path.join(directory, name) +'_split_randomized_rows_' + str(
                i) + '.tsv', header=False, index=False, sep='\t')
        if colnames is None:
            colnames = ['C' + str(x) for x in range(d.shape[1])]
        else:
            mycolnames = np.concatenate([np.asarray(['rownames']), colnames.values[start:(start+d.shape[1]-1)]])
            start= start+d.shape[1]-1
        d.columns = mycolnames
        d.to_csv(os.path.join(directory, name) +'_split_randomized_rows_colnames_' + str(
                i) + '.tsv', header=True, index=False, sep='\t')

def split_gwas():
    # could be made more realistic by adding actual snp names
    dirname = '/home/anne/Documents/featurecloud/dev/pca-tool/fed-pca-client/test_datasets/chr1'
    os.makedirs(dirname, exist_ok=True)
    data = pd.read_csv('/home/anne/Documents/featurecloud/pca/vertical-pca/data/1000g/raw/chr1/chr1.thin.traw.scaled', sep='\t')

    data = data.values
    data = data.astype('float')
    split_and_write(data, dirname, 'chr1')


if __name__ == '__main__':
    #split_mnist()
    #split_compass()
    split_gwas()





