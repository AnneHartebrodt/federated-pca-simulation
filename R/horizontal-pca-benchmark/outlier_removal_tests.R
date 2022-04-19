require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

robust.outlier<-function(x, m) which( (abs(x - median(x)) / mad(x)) > m )

read_and_pca<-function(type){
  data<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.tsv'))
  names<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.names.tsv'),
               header=F, sep='\t')
  
  names[,c('type', 'tss', 'id') := tstrsplit(V1, '-')]
  pca <- prcomp(data, scale. = F, center = F)
  
  x <- as.data.table(pca$x)[, 1:10]
  x$id<-names$id
  x$tss<-names$tss
  x$name<-names$V1
  return(x)
}

outdir <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/leave1out/sample/'

data_lo <- read_cancer_types(outdir = outdir)
data_lo[,c('type', 'tss', 'id') := tstrsplit(algorithm, '-')]

data_lo<- as.data.table(data_lo)
outlier_summarised<-data_lo[value>15 & name %in% paste0('', 1:10), .N, by=c('algorithm', 'dataset')]
ol<-outlier_summarised[N>2]