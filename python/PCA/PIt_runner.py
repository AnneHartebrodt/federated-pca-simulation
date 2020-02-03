from import_export.import_data import CustomDataImporter
from power_it import DistributedPowerIteration
import numpy as np
import scipy as sc
import os.path as path
import os as os

class PowerIterationRunner():
    def __init__(self, file = '/tmp/'):
        self.powerit = DistributedPowerIteration(file = file)
        self.importer = CustomDataImporter()

    def run_standalone(self, data, outfile=None, dims=1000, header=None, rownames=None, center=True, scale_var=True,scale01=False, scale_unit=True, transpose=False, sep='\t', filename='/pca', drop_samples=[], log=True, exp_var=0.5, epsilon=1, delta=1, noise=False):
        '''
        This function performs a regular principal component analysis and saves the result to files containing
        the projection the
        :param datafile: Unscaled datafile
        :param ddpca:
        :param outfile: path and name of the output file without extension
        :param dims: Number of dimensions to return (#eigenvectors and corresponding eigenvalues)
        :param seed: random seed
        :param nr_samples: #variables to select if not all the data columns are to be used for the pca
        :param header: row number which contains the header/ number of header rows
        :param rownames: column number which contains the rownames/sample ids
        :return: projection, eigenvectors and eigenvalues
        '''

        # if data is a string, it is the filename and has to be imported
        if isinstance(data, str):
            data, sample_ids, variable_names = self.importer.data_import(data, header=header, rownames=rownames,outfile=outfile, sep=sep, transpose=transpose)
        # remove outlier samples, previously identified
        if len(drop_samples) != 0:
            data = np.delete(data, drop_samples, 0)
        # drop columns with 0 variance or add some pseudocounts with low variance
        data, varn = self.importer.drop0Columns(data, None, noise=True, drop=False)
        # log(N+1) the data
        if log:
            data = self.importer.log_transform(data)
        # scale data
        data = self.importer.scale_data(data, center, scale_var, scale01, scale_unit)
        # run pca
        if noise:
            params = self.powerit.determine_parameters(cov, epsilon, delta, n)
            pca, W, s, noise = self.powerit.power_method(cov, params['sigma'], params['L'], params['p'], noise)
        else:
            noise_variance = 1
            L = 15
            p = 10
            pca, W, s, noise = self.powerit.power_method(cov, noise_variance, L, p, noise)
        if outfile is not None and filename is not None:
            self.save_PCA(pca, W, s, outfile + filename)
        # return the projected datapoints, the eigenvectors and the eigenvalues
        return pca, W, s, params


    def generate_result_str(self, inputfile, epsilon, delta, n, d, q, p, r, coherence, L, noise_variance, e, noise_vals, assumed_noise):
        res = path.basename(path.dirname(inputfile)) + '\t' + str(epsilon)+ '\t' + str(delta) + '\t' + str(n)+'\t' + str(d)+'\t' + str(q)+'\t' + str(p)\
              +'\t'+str(r)+'\t'+ str(coherence)+ '\t'+ str(L)+'\t' + str(noise_variance)+ '\t'+ str(e)+'\t'+str(assumed_noise)+'\t'
        header = 'study.id\tepsilon\tdelta\tnr.samples\tnr.dimension\tnr.itermediate\tit.rank\tr\tcoherence\tnr.iterations\tnoise.variance\terror.bound\tassumed.noise\t'
        i = 1
        for v in noise_vals:
            res = res +str(v) + '\t'
            header = header+'n1.'+str(i)+'\t'
            i= i+1
        return (res, header)

    def write_summary(self, res, header, outfile):
        try:
            os.makedirs(path.dirname(outfile))
        except OSError:
            print(path.dirname(outfile) + ' could not be created')
        else:
            print(path.dirname(outfile) + ' was created')

        if not path.exists(outfile):
            with open(outfile, 'w') as handle:
                handle.write(header + '\n')

        with open(outfile, 'a+') as handle:
            handle.write(res + '\n')


if __name__ == '__main__':
    # parser = ap.ArgumentParser(description='Eigengap calculation')
    # parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    # parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    # parser.add_argument('-e', metavar='epsilon', type=str, help='Epsilon, the privacy parameter')
    # parser.add_argument('-d', metavar='delta', type=str, help='Delta the privacy parameter')
    #
    # args = parser.parse_args()
    #
    # inputfile = args.f
    # outfile = args.o
    # epsilon = args.e
    # delta = args.d


    inputfile = '/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_only.tsv'
    outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/summary.txt'
    epsilon = 1
    delta = 0.1
    sep = '\t'


    DPIT = DistributedPowerIteration()

    # Scale data an calculate eigengap
    data, sample_id, var_names = DPIT.cd.data_import(inputfile, header=0, rownames=0, sep= sep)
    data_scaled = DPIT.cd.scale_data(data, center=True, scale_var=True, scale_unit=False, scale01=False)
    n = data.shape[0]
    cov = np.cov(data_scaled.T)

    pit = PowerIterationRunner()
    pca,W, E, params = pit.run_standalone(cov)
    print('finished')
    print(E[0:20])

    # res, header = DPIT.generate_result_str(inputfile, epsilon, delta, n, d, q1, p, r, coher, L, noise_variance, bound,noise, assumed_noise)
    #
    # DPIT.write_summary(res, header, outfile)
    #
    # pd.DataFrame(eigenvals).to_csv(path.dirname(inputfile) + '/eigenvalues_noisy_pit.tsv', sep='\t')
    # pd.DataFrame(eigenvectors).to_csv(path.dirname(inputfile) + '/eigenvectors_noisy_pit.tsv', sep='\t')
    #
    # #
    # eigenvectors, noise = DPIT.power_method(cov, noise_variance, int(L), q1, False)
    # eigenvals = []
    # for v in range(eigenvectors.shape[1]):
    #     eigenvals.append(DPIT.eigenvalue(cov, eigenvectors[:, v]))
    #
    # pd.DataFrame(eigenvals).to_csv(path.dirname(inputfile) + '/eigenvalues_pit.tsv', sep='\t')
    # pd.DataFrame(eigenvectors).to_csv(path.dirname(inputfile) + '/eigenvectors_pit.tsv', sep='\t')