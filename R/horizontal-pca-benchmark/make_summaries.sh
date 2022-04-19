#!/bin/bash

## Run in the root of the output directory


colname='angle'
outfile1='angles.u.tsv'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s 'angles.u' -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o 'angles.u.summary.tsv' -c $colname


outfile1='angles.v.tsv'
outfile2='angles.v.summary.tsv'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s 'angles.v' -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='conv.tsv'
outfile2='conv.summary.tsv' 
colname='conv'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s $colname -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname

outfile1='eigenvector_convergence.tsv'
outfile2='eigenvector_convergence.summary.tsv' 
colname='eigenvector_convergence'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s $colname -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='cor.tsv'
outfile2='cor.summary.tsv' 
colname='cor'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s $colname -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname


outfile1='eigenval.tsv'
outfile2='eigenval.summary.tsv' 
colname='eigenval'
Rscript ../../../federated_dp_pca/R/read_data.R -b . -s $colname -c $colname -o $outfile1
Rscript ../../../federated_dp_pca/R/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname
