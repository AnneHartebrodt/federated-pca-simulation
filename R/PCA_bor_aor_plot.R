if(!require(optparse)){
  install.packages('optparse')
}
require(data.table)
require(ggplot2)
require(ggbiplot)
require(bigutilsr)


option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type='character',
              help="coding genes output file"),
  make_option(c("-o", "--output.folder"), action="store", default=NA, type='character',
              help="The output file"),
  make_option(c("-k", "--number.clusters"), action="store", default=NA, type='numeric',
              help="The number of clusters for cluster detection")
)
opt = parse_args(OptionParser(option_list=option_list))

k = opt$k
pca.dir<-file.path(opt$o, 'pca')
dir.create(pca.dir)
lof.dir<-file.path(opt$o, 'lof')
dir.create(lof.dir)

#https://www.r-bloggers.com/detecting-outlier-samples-in-pca/
robust.outlier<-function(x) which( (abs(x - median(x)) / mad(x)) > 6 )
explained.variance<-function(x) cumsum(pca$sdev^2/sum(pca$sdev^2))
which.perc<-function(x, perc) min(which(x>=perc))

#Read and log transform data
#opt$f<-'/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_only.tsv'
data<-fread(file = opt$f)
var0<-which(apply(data,2, function(x) var(x)!=0))
data<-data[, unique(c(which(colSums(data)!=0), var0)), with=F]
data<-log2(data+1)

#run first PCA
pca<-prcomp(data, center = T, scale. = T)

for ( i in 2:10){
  ou<-unique(as.numeric(unlist(apply(pca$x[,1:i], 2, function(a) robust.outlier(a)))))
  lof = LOF(pca$x[, 1:i], seq_k = c(5,10,20))
  lof.hist<-ggplot(data.table(lof), aes(x = lof))+geom_histogram()
  ggsave(lof.hist, filename = file.path(lof.dir, paste0('hist_lof_', i, '.pdf')))
  lof<-which(lof>1)
  line<-paste(c(basename(dirname(opt$f)), length(ou), ou), collapse  = '\t')
  write(line,file=file.path(opt$o, 'outliers_mad.txt'),append=TRUE)
  line<-paste(c(basename(dirname(opt$f)), length(lof), lof), collapse  = '\t')
  write(line,file=file.path(opt$o, 'outliers_lof.txt'),append=TRUE)
}

ou<-unique(as.numeric(unlist(apply(pca$x[,1:3], 2, function(x) robust.outlier(x)))))
lof = LOF(pca$x[, 1:3], seq_k = c(5,10,20))
lof<-which(lof>1)



if(length(ou)>0){
groups <- rep('regular', nrow(data))
groups[ou]<-'mad'
groups[lof]<-'lof'
groups[intersect(ou, lof)]<-'all'
}else{
  groups <- rep('regular', nrow(data))
}

p<-ggbiplot(pca, var.axes = F, labels = 1:nrow(data), groups = as.factor(groups))
file.name<-paste0(basename(dirname(opt$f)), '_1_2.pdf')
ggsave(p, filename = file.path(pca.dir, file.name))

p<-ggbiplot(pca, var.axes = F, labels = 1:nrow(data), groups = as.factor(groups), choice = 2:3)
file.name<-paste0(basename(dirname(opt$f)), '_2_3.pdf')
ggsave(p, filename = file.path(pca.dir, file.name))

p<-ggbiplot(pca, var.axes = F, labels = 1:nrow(data), groups = as.factor(groups), choice = c(1,3))
file.name<-paste0(basename(dirname(opt$f)), '1_3.pdf')
ggsave(p, filename = file.path(pca.dir, file.name))

#Remove outlier and rerun.
if(length(unique(c(ou, lof)))!=0){
outlier.free<-data[-unique(c(ou, lof)), with=TRUE]
lab<-1:nrow(data)
lab<-lab[-unique(c(ou, lof))]
}else{
  outlier.free<-data
  lab<-1:nrow(data)
}
var0<-which(apply(outlier.free,2, function(x) var(x)!=0))
v<-unique(c(which(colSums(outlier.free)!=0), var0))
outlier.free<-outlier.free[, v , with=FALSE] 

pca.outlier.free<-prcomp(outlier.free, scale. = T, center = T)

p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab)
#p.out

p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab)
file.name.out<-paste0(basename(dirname(opt$f)),'_outlier_free_1_2', '.pdf')
ggsave(p.out, filename = file.path(pca.dir, file.name.out))

p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab, choices = c(1,3))
file.name.out<-paste0(basename(dirname(opt$f)),'_outlier_free_1_3', '.pdf')
ggsave(p.out, filename = file.path(pca.dir, file.name.out))

p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab, choices = c(2,3))
file.name.out<-paste0(basename(dirname(opt$f)),'_outlier_free_2_3', '.pdf')
ggsave(p.out, filename = file.path(pca.dir, file.name.out))


line <- paste0(sapply(seq(0.1, 1, by=0.1), function(x) which.perc(explained.variance(pca), x)), collapse = '\t')
write(line,file=file.path(opt$o, 'var_explained_bor.txt'),append=TRUE)
line<-paste0(sapply(seq(0.1, 1, by=0.1), function(x) which.perc(explained.variance(pca.outlier.free), x)), collapse = '\t')
write(line,file=file.path(opt$o, 'var_explained_aor.txt'),append=TRUE)


       
