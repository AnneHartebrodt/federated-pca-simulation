require(data.table)
require(ggplot2)
require(tidyr)

setwd('~/Documents/featurecloud/gwas/chr10//')



files <- c('fed_vs_plink.tsv','svd_vs_plink.tsv', 'qr_vs_plink.tsv', 'scipy_vs_fed.tsv',  'scipy_vs_svd.tsv', 'scipy_vs_qr.tsv')

data <- list()
for (file in files) {
  if (file.exists(file) & file.size(file)!=0){
  dat <- fread(file)
  dat$file <- gsub('.tsv', '', file)
  data[[file]]<-dat
  }
}
data<-rbindlist(data)

sub.rec.err<-data[,c(1:10, 21)]
sub.rec.err<- sub.rec.err %>% pivot_longer(-file)
sub.rec.err$name<-sapply(sub.rec.err$name, function(x) as.numeric(gsub('V', '', x)))
ggplot(sub.rec.err, aes(as.factor(name), value, fill=file))+
  geom_boxplot()+
  ggtitle('Angles between learning eigenvectors')+
  xlab('#eigenvector')+ylab('angle [degree]')+
  facet_wrap('file')


angles<-data[, c(11:20, 21)]
angles <- angles %>% pivot_longer(-file)
angles$name <- sapply(angles$name, function(x) as.numeric(gsub('V', '', x))-10)
ggplot(angles, aes(as.factor(name), value, fill=file))+geom_boxplot()+
  ggtitle('Mean of explained variances')+
  xlab('Eigenvalue rank')+ ylab('MEV')+
  facet_wrap('file')


eigenvalues <- list()
evals<- c('eigenvalues_fed.tsv', 'eigenvalues_qr.tsv', 'eigenvalues_svd.tsv')
for (ev in evals){
  edat <- fread(ev)
  edat$file <- ev
  eigenvalues[[ev]]<-edat
}
eigenvalues<-rbindlist(eigenvalues)

eigenvalues<- eigenvalues %>% pivot_longer(-file)
eigenvalues$name<-sapply(eigenvalues$name, function(x) as.numeric(gsub('V', '', x))-30)
ggplot(eigenvalues, aes(as.factor(name), value, fill=file))+
  geom_boxplot()+
  ggtitle('Eigenvalues')+
  facet_wrap('file')


mevs <- list()
evals<- c('sub_rec_err_fed.tsv', 'sub_rec_err_qr.tsv', 'sub_rec_err_svd.tsv')
for (ev in evals){
  edat <- fread(ev)
  edat$file <- ev
  mevs[[ev]]<-edat
}
mevs<-rbindlist(mevs)

mevs<- mevs %>% pivot_longer(-file)
mevs$name<-sapply(mevs$name, function(x) as.numeric(gsub('V', '', x)))

ggplot(mevs, aes(as.factor(name), value))+geom_boxplot()+
  ggtitle('Subspace reconstruction error')+
  xlab('# of eigenvectors used')+
  ylab('Error')+
  facet_wrap('file')
