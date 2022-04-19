basedir=/home/anne/Documents/featurecloud/pca
scriptpath=/home/anne/Documents/featurecloud/pca/federated_dp_pca/
export PYTHONPATH=$PYTHONPATH:$scriptpath

outfile=$basedir'/approximative-vertical/results'

echo "making summaries"
cd $outfile
echo $(pwd)


colname='angle'
outfile1='angles.u.tsv'
outfile2='angles.u.summary.tsv'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s 'angles.u' -c $colname -o $outfile1 -d $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2 -c $colname
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data_with_dummy.R -f $outfile1 -o $outfile2  -c $colname

outfile1='angles.v.tsv'
outfile2='angles.v.summary.tsv'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s 'angles.v' -c $colname -o $outfile1 -d $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='conv.tsv'
outfile2='conv.summary.tsv'
colname='conv'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname

outfile1='eigenvector_convergence.tsv'
outfile2='eigenvector_convergence.summary.tsv'
colname='eigenvector_convergence'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s $colname -c $colname -o $outfile1 $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='cor.tsv'
outfile2='cor.summary.tsv'
colname='cor'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='eigenval.tsv'
outfile2='eigenval.summary.tsv'
colname='eigenval'
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/read_data_2.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $scriptpath -m
Rscript $scriptpath/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


#bash $scriptpath/misc_scipts/make_summaries.sh $scriptpath

