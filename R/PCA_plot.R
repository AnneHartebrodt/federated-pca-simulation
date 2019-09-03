require(data.table)
require(ggplot2)
require(RColorBrewer)
require(viridisLite)
require(colorspace)

p1 = brewer.pal(10, 'Paired')
p2 = lighten(p1, 0.5)
p3 = darken(p1, 0.5)
pal = c(p1,p2,p3)


orig <- '/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/data_sample.tsv'
orig.data<- fread(orig, header = T, sep = '\t')
orig.data<-t(orig.data)
pcComp<- prcomp(orig.data, center = T, scale. = T)
regularPCA<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/standalonePCA.tsv')
regularPCA<-regularPCA[,0:2]
multiPCA1<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/multi_no_noise_1.tsv')
multiPCA2<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/multi_no_noise_2.tsv')
multiPCA3<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/multi_no_noise_3.tsv')
multiPCA4<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/multi_no_noise_4.tsv')
multiPCA5<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/multi_no_noise_5.tsv')

g1 <- ggplot(as.data.table(pcComp$x), aes(PC1,-PC2, col="#FDBF6F"))+geom_point()+
  ggtitle('PCA results R, python and python multisite')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))+
  theme(legend.key = element_blank())
g1

g1 <- g1 + geom_point(data = regularPCA, aes(x=V1, y = V2, col = "#A6CEE3"))
g1

g1<-g1 + geom_point(data = multiPCA1, aes(x=V1, y = V2, col = "#466E00"))
g1

g1<-g1 + geom_point(data = multiPCA2, aes(x=V1, y = V2, col = "#3C185E"))
g1
g1<-g1 + geom_point(data = multiPCA3, aes(x=V1, y = V2, col = "#815802"))
g1
g1<-g1 + geom_point(data = multiPCA4, aes(x=V1, y = V2, col = "#CAB2D6"))
g1
g1<-g1 + geom_point(data = multiPCA5, aes(x=V1, y = V2, col = "#C8F6A1"))
g1

pal

g1 <- ggplot(as.data.table(multiPCA2), aes(V1,V2, col="#466E00"))+geom_point()+
  ggtitle('PCA results R, python and python multisite')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))
g1




orig <- datasets::iris
#orig.data<- fread(orig, header = T, sep = '\t')
#orig.data<-t(orig.data)
#pcComp<- prcomp(orig[,1:4], center = T, scale. = T)
#orig <- '/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/data_sample.tsv'
#orig.data<- fread(orig, header = T, sep = '\t')
#orig<-t(orig.data)
pcComp<- prcomp(orig[,1:4], center = T, scale. = T)
orig<-as.data.frame(orig[1:4])
orig2 <-orig



orig2$max<-apply(orig2, c(1), function(x) max(x))
orig2$mim<-apply(orig2, c(1), function(x) min(x))

for (irow in 1:nrow(orig2)){
  for (jcol in 1:ncol(orig2[,1:4])){
    orig2[irow, jcol]<-orig2[irow, jcol]-orig2$mim[irow]/(orig2$max[irow]-orig2$mim[irow])
  }
}

pca2<-prcomp(orig2[,1:4], center = T, scale. = T)

dd= as.factor(orig$Species)
ggplot(as.data.table(pca2$x), aes(PC1, PC2, col='1'))+geom_point()+geom_point(data=as.data.table(pcComp$x), aes(PC1, -PC2, col = '2'))



orig <- datasets::iris
orig.data<- fread(orig, header = T, sep = '\t')
#orig.data<-t(orig.data)
#pcComp<- prcomp(orig[,1:4], center = T, scale. = T)
orig <- '/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/data_sample.tsv'
orig<- fread(orig, header = T, sep = '\t')
orig<-t(orig)
pcComp<- prcomp(orig, center = T, scale. = T)
orig<-as.data.frame(orig)
orig2 <-orig



orig2$max<-apply(orig2, c(1), function(x) max(x))
orig2$mim<-apply(orig2, c(1), function(x) min(x))

for (irow in 1:nrow(orig2)){
  for (jcol in 1:ncol(orig2[,1:1999])){
    orig2[irow, jcol]<-orig2[irow, jcol]-orig2$mim[irow]/(orig2$max[irow]-orig2$mim[irow])
  }
}

pca2<-prcomp(orig2[,1:1999], center = T, scale. = T)

dd= as.factor(orig$Species)
ggplot(as.data.table(pca2$x), aes(PC1, PC2, col='1'))+geom_point()+geom_point(data=as.data.table(pcComp$x), aes(PC1, -PC2, col = '2'))
