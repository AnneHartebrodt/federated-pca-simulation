basepath=$1
scriptpath=$2
export PYTHONPATH=$PYTHONPATH:$scriptpath

outfile=$basepath'/results-new-tests/'

mkdir -p $basepath/figures
Rscript $scriptpath/R/vertical-pca-benchmark/visualisation/iterations_vs_convergence.R -f '$outfile/1000g/chr1.tsv,$outfile/1000g/chr2.tsv,$outfile/mnist/long.dummy.angles.u.summary.tsv' \
 -o $basepath/figures/angles_all.pdf -n 'Chromosome 1,Chromosome 2,MNIST'
