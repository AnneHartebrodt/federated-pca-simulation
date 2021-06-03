from python.PCA.vertical.vertical_pca_benchmark import *
import argparse as ap


def run_scalability_test(data, dataset_name,outdir):
    for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    #for p in [0.6, 0.7, 0.8, 0.9]:
        od = op.join(outdir, str(p))
        os.makedirs(od, exist_ok=True)
        # repeat 10 times
        for i in range(10):
            # select arbitrary subset of samples
            data_list, choices = sh.partition_data_vertically(data, splits=2, randomize=True, perc=[p, 1-p], equal=False)
            # run all algorithms once
            the_epic_loop(data=data_list[0], dataset_name=dataset_name, maxit=500, nr_repeats=1, k=10, splits=[3],
                          outdir=od, ortho_freq=100)

if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Split datasets and run "federated PCA"')
    parser.add_argument('-f', metavar='file', type=str, help='mnist input dir')
    parser.add_argument('-o', metavar='file', type=str, help='output dir')
    args = parser.parse_args()


    # read data
    input_dir = args.f
    #data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)
    data = si.scale_center_data_columnwise(data, center=True, scale_variance=False)
    # Images=samples are vertically distributed
    data = data.T
    dataset_name = 'mnist'
    #outdir = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/scalability2'
    outdir = args.o
    run_scalability_test(data, dataset_name, outdir=outdir)


