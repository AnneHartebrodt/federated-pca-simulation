require(data.table)
require(ggplot2)

standalone<-fread('/home/anne/Documents/featurecloud/data/TCGA/htseq/simulation/output.pca')

ggplot(standalone, aes(x=V1, y=V2))+geom_point()+geom_text(aes(label=rownames(standalone)),hjust=0, vjust=0)+ggtitle('Standalone all')

standalone.no.out<-standalone[-c(863,952)]

ggplot(standalone.no.out, aes(x=V1, y=V2))+geom_point()+
  geom_text(aes(label=rownames(standalone)[-c(863,952)]),hjust=0, vjust=0)+
  ggtitle('Standalone without outliers')


data <- fread('/home/anne/Documents/featurecloud/data/TCGA/htseq/simulation/eigenvalues/projections.txt')

cn<-1:(ncol(data)-4)
cn<-paste0('PC', cn)


nrdat<-nrow(unique(data[,c(ncol(data)-3,ncol(data)-2,ncol(data)-1,ncol(data)), with=F]))
indx<-nrow(data)/nrdat

data$index <- rep(1:indx,nrdat)
ggplot(data, aes(x=V1, y=V2))+geom_point()+geom_text(aes(label=index),hjust=0, vjust=0)

data.no.out<-data[!index %in% c(863,952)]
ggplot(data.no.out[V10==100 & V9==0 & V8==1], aes(x=V1, y=V2))+geom_point()+geom_text(aes(label=index),hjust=0, vjust=0)

