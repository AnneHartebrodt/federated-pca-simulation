# Scripts for benchmark of vertically-federated Principal Component Analysis

This is the workflow followed to benchmark federated vertically partitioned PCA. 

Careful: GWAS data can be large and some steps might not run (or at least not very fast) on a standard desktop computer. Since this is a simulation environment, the tests are repeated several times and therefor it will run a while.


## Create an environment
```
basepath=$pwd/vertical-pca
mkdir -p $basepath
cd $basedir
git clone https://gitlab.com/hartebrodt/federated_dp_pca.git --recursive
scriptpath=$pwd/federated_dp_pca
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

outfile=$basepath'/results/mnist'


mkdir -p $outfile
python3 $scriptpath/python/PCA/vertical/vertical_pca_benchmark.py -f $basepath'/data/mnist/raw' \
--filetype 'mnist' --center -o $outfile -r 20 -k 10 \
 -i 2000 -s 5 --vert --ortho_freq 100

echo "making summaries"
cd $outfile
echo $(pwd)
bash $scriptpath/misc_scipts/make_summaries.sh $scriptpath
```
## run GWAS tests
### get plink1 and 2
```
mkdir -p $basepath/software
cd $basepath/software
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20210608.zip
unzip plink_linux_x86_64_20210606.zip
unzip plink2_linux_x86_64_20210608.zip
plink1path=$basepath/software/plink
plink2path=$basepath/software/plink2
cd ..


### Get the data
```
bash $scriptpath/misc_scipts/1000g_download.sh $pwd $plink1path $plink1path
```

### Run GWAS tests
```
bash $scriptpath/misc_scipts/run_tests_gwas.sh
```


## visualize the results
```
mkdir -p $basepath/figures
Rscript $scriptpath/R/vertical-pca-benchmark/visualisation/iterations_vs_convergence.R -f '$outfile/1000g/chr1.tsv,$outfile/1000g/chr2.tsv,$outfile/mnist/long.dummy.angles.u.summary.tsv' \
 -o $basepath/figures/angles_all.pdf -n 'Chromosome 1,Chromosome 2,MNIST'

# MNIST only
#mkdir -p $basepath/figures
#Rscript $scriptpath/R/vertical-pca-benchmark/visualisation/iterations_vs_convergence.R -f '$outfile/mnist/long.dummy.angles.u.summary.tsv' -o $basepath/figures/angles_all.pdf -n 'MNIST'

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
