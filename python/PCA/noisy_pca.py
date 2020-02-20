import argparse as ap
import numpy as np
import scipy as sc
import convenience as cv
import comparison as co
import os.path as path
import proxy_covariance_runner as runner
from import_export.import_data import CustomDataImporter
import outlier_removal as orm
import power_iteration as power

def compare_noisy_with_regular(datafile, outfile=None, dims=1000, header=0, rownames=0, center=True, scale_var=True, scale01=False, scale_unit=True,transpose=False, sep = '\t', reported_angles = 20, exp_var = 0.5, log=True):

    study_id = path.basename(path.dirname(datafile))

    importer = CustomDataImporter()
    data, varnames, sampleids = importer.data_import(datafile, header=header, rownames=rownames, outfile=outfile,
                                                          transpose=False, sep=sep)
    n = data.shape[0]
    dims = min(dims, n)

    pca, W1, E1 = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,
                                                 center=center, scale_var=scale_var, scale01=scale01,
                                                 scale_unit=scale_unit, transpose=transpose, sep=sep,
                                                 filename='/pca.before_outlier_removal', log=log, exp_var=exp_var)
    outliers = orm.outlier_removal_mad(pca)

    pca_of, W1_of, E1_of = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,
                                     center=center, scale_var=scale_var, scale01=scale01,
                                     scale_unit=scale_unit, transpose=transpose, sep=sep,
                                     filename='/pca.before_outlier_removal', log=log, exp_var=exp_var, drop_samples=outliers)

    epsilons = [0.01, 0.05, 0.1, 0.5, 1, 5, 10, 100]
    deltas = [0.01]

    for epsilon in epsilons:
        for delta in deltas:
            pca_noisy, W1_noisy, E1_noisy = runner.run_standalone(data, outfile, dims=dims, header=header, rownames=rownames,center=center, scale_var=scale_var, scale01=scale01,scale_unit=scale_unit, transpose=transpose, sep=sep,filename='/pca.before_outlier_removal', log=log,exp_var=exp_var, noise=True, epsilon=epsilon, delta=delta, drop_samples=outliers)

            angles = co.compute_angles(W1_of, W1_noisy, reported_angles=reported_angles)
            with open(path.join(outfile,'angles_regular_noisy.tsv'), 'a+') as handle:
                id = study_id+'\t'+str(epsilon)+'\t'+str(delta)
                handle.write(cv.collapse_array_to_string(angles, study_id=id))

if __name__=="__main__":
    print('run comparison script')

    # parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    # parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    # parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default = 100)
    # args = parser.parse_args()
    #
    # inputfile = args.f
    # outfile = args.o
    # sep = args.s
    # dims = args.d


    compare_noisy_with_regular(inputfile, outfile, dims = dims, scale_unit=True, sep = sep, reported_angles=20, exp_var =0.99, rownames=None, header = 0, log = False)