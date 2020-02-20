import argparse as ap
import pandas as pd
import import_export.easy_import as easy
import proxy_covariance as dpca


def run_standalone(data, outfile=None, dims=100, header=None, rownames=None, center=False, scale_var=False,scale01=False, scale_unit=False, transpose=False, sep='\t', filename=None, drop_samples=[],log=False, epsilon=1, delta=1, noise=False):
    """
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
    """

    # if data is a string, data needs to be read first, otherwise it is
    # assumed to be scaled and ready to use
    if isinstance(data, str):
        data = easy.easy_import(data, header=header, rownames=rownames, sep=sep, center=center, scale_var=scale_var,scale01=scale01, scale_unit=scale_unit, drop_samples=drop_samples, log=log, outfile=outfile,transpose=transpose)

    pca, W, s = dpca.standalone_pca(data, ndims=dims, noise=noise, epsilon=epsilon, delta=delta)
    if outfile is not None and filename is not None:
        save_PCA(pca, W, s, outfile + filename)
    # return the projected datapoints, the eigenvectors and the eigenvalues
    return pca, W, s


def save_PCA(pca, W, s, outfile):
    if pca is not None:
        pd.DataFrame(pca).to_csv(outfile + '.projection', sep='\t', header=None, index=False)
    if W is not None:
        pd.DataFrame(W).to_csv(outfile + '.eigenvectors', sep='\t', header=None, index=False)
    if s is not None:
        pd.DataFrame(s).to_csv(outfile + '.eigenvalues', sep='\t', header=None, index=False)


if __name__ == "__main__":
    print('run')

    parser = ap.ArgumentParser(description='Run distributed PCA simulation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-d', metavar='dimensions', type=int, help='number of principal components to return')
    parser.add_argument('-p', metavar='output directory', type=str, help='output directory for simulation study.')
    parser.add_argument('-k', metavar='number_hospitals', type=int, help='Number of simulated hospitals', default=5)
    parser.add_argument('-s', action='store_true',
                        help='If true the generated eigenspaces are saved (!a lot of large files!)', default=False)
    parser.add_argument('-r', metavar='Repeats', type=int, help='Number of repetitions of the sampling process',
                        default=5)
    parser.add_argument('-c', help='True of data has column headers',
                        default=False, action='store_true')
    parser.add_argument('-i', metavar='sampleids', type=int,
                        help='Dataframe column which contains the sample ids, 0 index,', default=-1)
    parser.add_argument('-e', metavar='epsilons', type=str, help='epsilons to simulate separated by a comma')
    parser.add_argument('-g', metavar='deltas', type=str, help='deltas to simulate separated by a comma')
    parser.add_argument('-v', action='store_true',
                        help='Scale variables to unit variance', default=True)
    parser.add_argument('-u', action='store_true',
                        help='Scale samples to unit norm', default=True)
    parser.add_argument('-z', action='store_true',
                        help='Scale variables between 0 and 1', default=False)
    parser.add_argument('-t', action='store_true',
                        help='Center variables by substracting the mean', default=True)

    parser.add_argument('-A', action='store_true',
                        help='Run standalone simulation', default=False)
    parser.add_argument('-B', action='store_true',
                        help='Run distributed simulation without noise', default=False)
    parser.add_argument('-C', action='store_true',
                        help='Run distributed simulation with noise', default=False)
    parser.add_argument('-D', action='store_true',
                        help='Run distributed, submit a list of files (comma sep), noise = F', default=False)
    parser.add_argument('-E', action='store_true',
                        help='Run distributed, submit a list of files (comma sep), noise = T', default=False)
    parser.add_argument('-V', action='store_true',
                        help='Save Covariance matrix before and after noise', default=False)
    args = parser.parse_args()


    def parse_array(value_str):
        values_str = value_str.split(',')
        values_int = []
        for v in values_str:
            values_int.append(float(v))
        return values_int


    if args.c:
        header = 0
    else:
        header = None

    if args.i == -1:
        rownames = None
    else:
        rownames = args.i

    try:
        epsilons = parse_array(args.e)
    except:
        print('Incompatible epsilon parameters')
        print('default to epsilon=0.1')
        epsilons = [0.1]

    try:
        deltas = parse_array(args.g)
    except:
        print('Incompatible delta parameters')
        print('default to delta=0.001')
        deltas = [0.001]

    print('file: ' + str(args.f))
    print('dimensions: ' + str(args.d))
    print('output directory: ' + str(args.p))
    print('number of hospitals: ' + str(args.k))
    print('save eigenvalues: ' + str(args.s))
    print('repeats: ' + str(args.r))
    print('column headers: ' + str(args.c))
    print('sample ids: ' + str(args.i))

    if args.A:
        run_standalone(args.f, outfile=args.p, dims=args.d, header=header, rownames=rownames,
                       center=args.t, scale_var=args.v, scale01=args.z, scale_unit=args.u,
                       transpose=False)