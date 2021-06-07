basepath=$1
scriptpath=$2
export PYTHONPATH=$PYTHONPATH:$scriptpath
echo $PYTHONPATH


outdir=$basepath'/results/scalability'
mkdir -p $outdir
indir=$basepath'/data/mnist/raw'
#python3 $scriptpath/python/PCA/vertical/scalability.py -f $indir -o $outdir

cd $outdir
echo $scriptpath/R/vertical-pca-benchmark/data_cleanup/scalability_aggregation.R
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/scalability_aggregation.R -b . -o 'transmission_cost.tsv' -d $scriptpath
echo 'making angle summaries'

colname='angle'
outfile1='angles.u.tsv'
outfile2='angles.u.summary.tsv'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data.R -b . -s 'angles.u' -c $colname -o $outfile1 -d $scriptpath -a
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2 -c $colname
