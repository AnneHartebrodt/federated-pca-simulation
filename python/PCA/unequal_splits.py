from master import Distributed_DP_PCA
from PCA_runner import SimulationRunner
from outlier_removal import OutlierRemoval
from import_export.import_data import CustomDataImporter
import os.path as path
import numpy as np
import convenience as cv
import comparison as co
import argparse as ap


class AngleRunner():
    def __init__(self, file = '/tmp/'):
        self.ddppca = Distributed_DP_PCA(file = file)
        self.importer = CustomDataImporter()
        self.simulation = SimulationRunner()
        self.outlier = OutlierRemoval()



    def unqeal_split(self, data, interval_end, scale_variance=True, center=True, scale01=False, scale_unit=False, ndims=1000, header=0, rownames=0,  scale_var=True,transpose=False, sep = '\t', exp_var = 0.5, mult_dims_ret =1):
        """
        This function simulates a multisite PCA with each site having
        varying number of samples.
        :param data:
        :param sites: number of sites to split the data into
        :return:
        """

        print('Shuffling data')
        np.random.shuffle(data)

        Ac = []
        vexx = []
        s = 0
        start = 0

        for i in range(len(interval_end)):
            s = s + 1
            end = int(interval_end[i])
            # slice matrix
            data_sub = data[start:end, :]
            
            print('Local PCA for outlier identification')
            pca, W1, E1 = self.simulation.run_standalone(data_sub, outfile, dims=ndims, header=header, rownames=rownames,center=center,scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.loc', drop_samples=[], log=True, exp_var=exp_var)

            print('Outlier removal')
            outliers = self.outlier.outlier_removal_mad(pca, 6, 3)
            if len(outliers) != 0:
                data_sub = np.delete(data_sub, outliers, 0)


            print('Drop empty variable, scale local data, compute covariance')
            data_sub, vn = self.importer.drop0Columns(data_sub,None, drop = False, noise=True)
            data_sub = self.importer.log_transform(data_sub)
            data_sub = self.importer.scale_data(data_sub, center=center, scale_var=scale_variance, scale01=scale01,
                                                scale_unit=scale_unit)
            print('calculating cov')
            noisy_cov = self.ddppca.compute_noisy_cov(data_sub, epsilon0=1, delta0=1, nrSamples=data.shape[0],
                                                      nrSites=len(interval_end), noise=False)
            start = int(interval_end[i])
            print('finished')

            #return the local PCA with the maximal number of dimensions required.
            m  = max(mult_dims_ret)
            A, vex = self.ddppca.perform_SVD(noisy_cov, ndims=ndims, mult_dims_returned=m, var_exp=exp_var)
            Ac.append(A)
            vexx.append(vex)

        Ws = []
        Xs = []
        vex = max(vexx)
        print('Aggregate local PCAs')
        for mm in mult_dims_ret:
            mm = int(np.ceil(vex*mm))
            W, X = self.ddppca.aggregate_partial_SVDs(Ac, ndim=ndims, intermediate_dims=mm)
            W = self.ddppca.normalize_eigenvectors(W)
            Ws.append(W)
            Xs.append(X)
        return Ws, Xs


    def run_and_compare_unequal(self, datafile, outfile=None, dims=1000, header=0, rownames=0, center=True, scale_var=True, scale01=False, scale_unit=False,transpose=False, sep = '\t', reported_angles = 20, exp_var = 0.5, mult_dims_ret = [1,2,1.5, 5]):

        # ID is the folder name
        study_id = path.basename(path.dirname(datafile))
        print('Standalone PCA before outlier removal')

        data, varnames, sampleids = self.importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile, transpose=False, sep=sep)
        n = data.shape[0]
        dims = min(dims, n)

        pca, W1, E1 = self.simulation.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames, center=center, scale_var=scale_var, scale01=scale01,scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.before_outlier_removal', log = True,exp_var=exp_var)


        print('Logging outliers')
        outliers = self.outlier.outlier_removal_mad(pca, 6, 3)
        print(outliers)
        with open(outfile + '/removed_outliers.tsv', 'a+') as handle:
            handle.write(cv.collapse_array_to_string(outliers, study_id))

        print('Standalone PCA after outlier removal')
        pca, W1, E1 = self.simulation.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames, center=center, scale_var=scale_var, scale01=scale01, scale_unit=scale_unit,transpose=transpose, sep=sep, filename='/pca.after_outlier_removal',drop_samples=outliers, log = True, exp_var=exp_var)

        with open(outfile + '/nr_vars_explain_aor.tsv', 'a+') as handle:
            for var in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
                handle.write(str(var)+'\t'+str(self.ddppca.variance_explained(E1, perc=var))+'\n')

        W1 = self.ddppca.normalize_eigenvectors(W1)

        interval_end = self.make_test_intervals(n)


        for ar in interval_end:
            for i in range(10): # run it more often the smaller the splits
                print('Current split')
                print(ar)
                Ws, Xs = self.unqeal_split(data, ar, ndims=dims, exp_var = exp_var, mult_dims_ret=mult_dims_ret)
                for w in range(len(Ws)) :
                    angles = self.compute_angles(W1, Ws[w], reported_angles=reported_angles)
                    with open(outfile + '/angles_unequal_splits'+str(mult_dims_ret[w])+'.tsv', 'a+') as handle:
                        handle.write(cv.collapse_array_to_string(angles, study_id=study_id))
                    with open(outfile + '/eigenvalues'+str(mult_dims_ret[w])+'.tsv', 'a+') as handle:
                        handle.write(cv.collapse_array_to_string(Xs[w][0:reported_angles], str(i)))

                    #write results
                meta = [len(ar)] + ar

                with open(outfile + '/meta_splits.tsv', 'a+') as handle:
                    handle.write(cv.collapse_array_to_string(meta,  str(i)))



    def make_test_intervals(self, n):
        # hardcoded for now
        # more extreme cases maybe later
        unequal_splits = list()
        unequal_splits.append([[1.0]])
        unequal_splits.append([[0.1, 0.9], [0.3, 0.7], [0.5, 0.5]])
        unequal_splits.append([[0.2, 0.2, 0.2, 0.2, 0.2], [0.1, 0.1, 0.2, 0.2, 0.4], [0.1, 0.1, 0.1, 0.1, 0.6],
                               [0.2375, 0.2375, 0.2375, 0.2375, 0.05]])

        interval_end = list()
        sum = 0
        for i in unequal_splits:
            for j in i:
                inter = list()
                for k in j:
                    sum = sum + k
                    inter.append(np.ceil(n * sum))
                sum = 0
                interval_end.append(inter)
        return interval_end

    def compute_angles(self, canonical, split, reported_angles=20):
        angles = list()
        for i in range(min(reported_angles, min(canonical.shape[1], split.shape[1]))):
            angle = co.angle(canonical[:, i], split[:, i])
            angles.append(angle)
        return(angles)

def parse_array(value_str):
    values_str = value_str.split(',')
    values_int = []
    for v in values_str:
        values_int.append(float(v))
    return values_int


if __name__=="__main__":
    print('run split script')

    #parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    #parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    #parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    #parser.add_argument('-v', metavar='explained_var', type=float, help='explained variance')
    #parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    #parser.add_argument('-m', metavar='mult_dims_ret', type=str, help='comma separated list of intermediate dimensions', default = 1)
    #parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default = 100)
    #args = parser.parse_args()

    #inputfile = args.f
    #outfile = args.o
    #exp_var = args.v
    #mult_dims_ret = args.m
    #sep = args.s
    #dims = args.d

    inputfile ='/home/anne/Downloads/mnist/train_flat.csv'
    outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/mnistt/'
    exp_var = 0.5
    sep = ','
    mult_dims_ret = '1,2,1.5,5'
    dims = 100


    sim = AngleRunner()

    mult_dims_ret = parse_array(mult_dims_ret)
    summaryfile = cv.make_eigenvector_path(outfile, path.basename(path.dirname(inputfile)))
    sim.run_and_compare_unequal(inputfile, summaryfile, dims = dims, scale_unit=False, sep = sep, reported_angles=20, exp_var =exp_var,
                                mult_dims_ret=mult_dims_ret, rownames = None, header = 0)

