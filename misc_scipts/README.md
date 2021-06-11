# Scripts for benchmark of vertically-federated Principal Component Analysis

## Create an environment
```
basepath=$pwd/vertical-pca
mkdir -p $basepath
cd $basedir
git clone https://gitlab.com/hartebrodt/federated_dp_pca.git --recursive
scriptpath=$pwd/federated
```

## Set up environment
```
cd federated_dp_pca
conda env create -f ../environment.yml
conda activate federated_dp_pca
cd ..
``

## Get the data
``
mkdir -p data/mnist/raw
cd data/mnist/raw
#wget http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz
#wget http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz
#wget http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz
#wget http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz
cd $basedir
```

## Run mnist test

```
export PYTHONPATH=$PYTHONPATH:$scriptpath

outfile=$basepath'/results-new-tests/mnist'


mkdir -p $outfile
python3 $scriptpath/python/PCA/vertical/vertical_pca_benchmark.py -f $basepath'/data/mnist/raw' \
--filetype 'mnist' --center -o $outfile -r 20 -k 10 \
 -i 2000 -s 5 --vert --ortho_freq 100

echo "making summaries"
cd $outfile
echo $(pwd)
bash $scriptpath/misc_scipts/make_summaries.sh $scriptpath
```

## visualize the results
```
mkdir -p $basepath/figures
Rscript $scriptpath/R/vertical-pca-benchmark/visualisation/iterations_vs_convergence.R -f '$outfile/1000g/chr1.tsv,$outfile/1000g/chr2.tsv,$outfile/mnist/long.dummy.angles.u.summary.tsv' \
 -o $basepath/figures/angles_all.pdf -n 'Chromosome 1,Chromosome 2,MNIST'

```

## Run scalability tests
```
outdir=$basepath'/results/scalability'
mkdir -p $outdir
indir=$basepath'/data/mnist/raw'
#python3 $scriptpath/python/PCA/vertical/scalability.py -f $indir -o $outdir

cd $outdir
echo $scriptpath/R/vertical-pca-benchmark/data_cleanup/scalability_aggregation.R
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/scalability_aggregation.R -b . -o 'transmission_cost.tsv' -d $scriptpath
echo 'making angle summaries'

mkdir -p $basepath/figures
Rscript $scriptpath/R/vertical-pca-benchmark/visualisation/scalability_vis.R -f  $basepath/results/scalability/summary.transmission_cost.tsv -o $basepath/figures/transmission_cost.pdf
```
