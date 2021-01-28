#!/bin/bash

## Run in the root of the output directory

## supply the basedir of the root scripts repo
if [ $# -eq 0 ]
  then
    basedir='/home/anne/Documents/featurecloud/pca/federated_dp_pca/'
  else
    echo 'basedir supplied'
    echo $1
    basedir=$1
fi



colname='angle'
outfile1='angles.u.tsv'
outfile2='angles.u.summary.tsv'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s 'angles.u' -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2 -c $colname


outfile1='angles.v.tsv'
outfile2='angles.v.summary.tsv'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s 'angles.v' -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='conv.tsv'
outfile2='conv.summary.tsv'
colname='conv'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname

outfile1='eigenvector_convergence.tsv'
outfile2='eigenvector_convergence.summary.tsv'
colname='eigenvector_convergence'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='cor.tsv'
outfile2='cor.summary.tsv'
colname='cor'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='eigenval.tsv'
outfile2='eigenval.summary.tsv' 
colname='eigenval'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname

outfile1='angles_precomp.tsv'
outfile2='angles_precomp.summary.tsv'
colname='angle'
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/read_data.R -b $(pwd) -s $colname -c $colname -o $outfile1 -d $basedir
Rscript $basedir/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname
