import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *
import math


def simulate_chen(data_list, uu, k=10, T=10,TT=5):
    '''

    Args:
        data_list: Data splits
        uu: canonical eigenvector
        k: subspace dimension
        T: NUmber of outer iterations
        TT: Number of inner iterations

    Returns:

    '''
    ###
    # Implementation of Chen et als
    # Distributed Estimation for Principal Component Analysis: a Gap-free Approach
    # Terrible convergence behavior
    ###
    u,s,v = lsa.svds(data_list[0], k=k)
    s = np.flip(s)
    eta = math.pow(data_list[0].shape[1]/data_list[0].shape[0], 0.25)
    s1 = np.sqrt(s[0])+ 1.5* eta
    print(s)
    Hklist = []
    for d in range(len(data_list)):
        Hk = s1 * np.identity(data_list[d].shape[1]) - np.dot(data_list[d].T, data_list[d]) / data_list[d].shape[0]
        Hklist.append(Hk)
    u,s,v = lsa.svds(Hklist[0])
    wk = np.flip(v.T, axis=0)
    Hinv = np.linalg.inv(Hklist[0])
    #Hinv = np.linalg.pinv(Hklist[0])
    wk = wk[:,0]
    # wkt1 = current
    wkt1 = wk
    conv = []
    for t in range(T):
        for t1 in range(TT):
            gklist = []
            for d in range(len(data_list)):
                gk = np.dot(Hklist[d], wkt1) - wk
                gklist.append(gk)
            # line 10
            g = (1/len(data_list))* np.nansum(gklist, axis=0)
            wkt1 = wkt1 - np.dot(Hinv, g)
        # line 12
        wkt1 = wkt1/np.linalg.norm(wkt1)
        wk = wkt1
        conv.append(co.angle(uu, wk))
    return wkt1, conv


def generate_data(number_features, number_samples, variable_means=None, variable_stds=None, seed=11):
    """
    Returns generated gaussian data either standard normal, if no means and standard deviations
    are given, or according to the given parameters.
    :param number_features: Number of variables to generate
    :param number_samples: Number of samples to generate
    :param variable_means: Optional vector of means for variables
    :param variable_stds: Optional vector of standard deviations for variables
    :param seed: random seed
    :return: numpy array of generated data containing samples (rows) x features (colums)
    """
    np.random.seed(seed)
    if variable_means is not None:
        assert  len(variable_means) == len(variable_stds), "Number of means and stds should be equal"
        assert len(variable_means) == number_features, "Number of variables should be eqaul to length of mean/std vectors"
    generated_data = []
    for nf in range(number_features):
        if variable_means is None:
            generated_data.append(np.random.standard_normal(number_samples))
        else:
            generated_data.append(np.random.normal(size=number_samples, loc=variable_means[nf], scale=variable_stds[nf]))
    generated_data = np.stack(generated_data, axis=1)
    return generated_data

if __name__ == '__main__':
    from python.PCA.horizontal.horizontal_pca_benchmark import read_presplit_data_folders, compute_canonical,scale_datasets

    import matplotlib.pyplot as pyp
    import seaborn as sns
    # generate data according to Chen et al.
    convergence_list  = []
    for i in range(5):
        data = np.random.normal(1, 0, 200*50*50).reshape( 200*50,50)
        # Eigengaps
        diag = [7,5,3,1,1]+[1]*45
        data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
        data_list, choices = sh.partition_data_horizontally(data, splits=200, randomize=False)
        # Orthonormalise matrices
        data_list = [la.qr(d, mode='economic')[0] for d in data_list]
        data_list = [np.dot(np.dot(d, np.diag(diag)), d.T) for d in data_list]
        data = np.concatenate(data_list, axis=0)
         #Compute canonical solution
        uu, ss, vv = lsa.svds(data, k=10)
        vv = np.flip(vv.T, axis=1)
        u, conv = simulate_chen(data_list, vv[:, 0], k=10,T=300, TT=5)
        convergence_list.append(conv)
        # Does converge on small simulated data
    convergence_list = np.stack(convergence_list, axis=0)
    mean = np.nanmean(convergence_list, axis=0)
    sd = np.nanstd(convergence_list, axis=0)
    sns.lineplot(x=range(len(mean)),y=mean, ci=sd)
    pyp.show()

    # Same test
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    data_list, choices = sh.partition_data_horizontally(data, splits=200, randomize=False)

    uu, ss, vv = lsa.svds(data, k=10)
    vv = np.flip(vv.T, axis=1)
    u, conv = simulate_chen(data_list, vv[:, 0], k=10,T=300, TT=5)

    # Does not converge within reasonable number of iterations
    pyp.plot(conv)
    pyp.show()

