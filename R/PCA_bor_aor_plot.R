if(!require(optparse)){
  install.packages('optparse')
}
require(data.table)
require(ggplot2)
require(ggbiplot)


option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type='character',
              help="coding genes output file"),
  make_option(c("-o", "--output.folder"), action="store", default=NA, type='character',
              help="The output file")
)
opt = parse_args(OptionParser(option_list=option_list))


robust.outlier<-function(x) which( (abs(x - median(x)) / mad(x)) > 6 )
explained.variance<-function(x) cumsum(pca$sdev^2/sum(pca$sdev^2))
which.perc<-function(x, perc) min(which(x>=perc))

#Read and log transform data
#opt$f<-'/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_trunc.tsv'
data<-fread(file = opt$f)
data<-data[, which(colSums(data)!=0), with=F]
data<-log2(data+1)

#run first PCA
pca<-prcomp(data, center = T, scale. = T)
ou<-robust.outlier(pca$x[,1:3])
p<-ggbiplot(pca, var.axes = F, labels = 1:nrow(data))
file.name<-paste0(basename(dirname(opt$f)), '.pdf')
ggsave(p, filename = file.path(opt$o, file.name))

#Remove outlier and rerun.
outlier.free<-data[-ou]
outlier.free<-outlier.free[, which(colSums(outlier.free)!=0), with=F] 
pca.outlier.free<-prcomp(outlier.free, scale. = T, center = T)
lab<-1:nrow(data)
p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab[-ou])
file.name.out<-paste0(basename(dirname(opt$f)),'_outlier_free', '.pdf')
ggsave(p.out, filename = file.path(opt$o, file.name.out))
