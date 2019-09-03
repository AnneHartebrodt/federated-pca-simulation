require(data.table)
require(ggplot2)
require(RColorBrewer)
require(viridisLite)
require(colorspace)
library(ggbiplot)
require(Matrix)
require(ggforce)

p1 = brewer.pal(10, 'Paired')
p2 = lighten(p1, 0.5)
p3 = darken(p1, 0.5)
pal = c(p1,p2,p3)
simu= fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/radition_simulation_results.tsv', sep = '\t', header = F)
#colnames(simu)= c('x', 'y', 'a', 'b', 'run', 'epsilon', 'delta', 'name')
#[run %in% c(5,6,7,8,9)& epsilon %in% c(0.001, 0.01, 0.05, 0.1, 0.5)]
g <- ggplot(simu[run %in% c(5,6,7,8,9)& epsilon %in% c(0.001, 0.01, 0.05, 0.1, 0.5)], aes(x, y, col=name))+geom_point()+facet_grid(c('epsilon', 'run'))+
  scale_color_manual(values = pal)+
  ggtitle('PCA Iris dataset comparing 5 runs and 5 epsilons')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))
plot(g)

simu= fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/radition_simulation_results.tsv', sep = '\t', header = F)
#colnames(simu)= c('x', 'y', 'a', 'b', 'run', 'epsilon', 'delta', 'name')
#[run %in% c(5,6,7,8,9)& epsilon %in% c(0.001, 0.01, 0.05, 0.1, 0.5)]
g <- ggplot(simu[80002:10001], aes(V1, V2))+geom_point()+
  scale_color_manual(values = pal)+
  ggtitle('PCA Iris dataset comparing 5 runs and 5 epsilons')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))
plot(g)

regularPCA<-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/radition_PCA_centered.tsv')
g1 <- ggplot(regularPCA, aes(V1,V2))+geom_point()+
  scale_color_manual(values = pal)+
  ggtitle('PCA Iris dataset comparing 5 runs and 5 epsilons')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))




orig<-fread('/home/anne/Documents/featurecloud/data/radiation/data.tsv')
orig <- as.data.table(transpose(orig))

pca <- svd(orig)
w <-pca$v

scalar1 <- function(x) {x / sqrt(sum(x^2))}

rot<-as.matrix(w[,1:2])
rot[,1]<--rot[,1]
rot[,2]<-rot[,2]
rot<-as.matrix(rot)
a <- scalar1(rot[,1])
plo<-as.matrix(orig) %*% rot

g3 <- ggplot(data.table(plo), aes(V1,V2))+geom_point()+
  scale_color_manual(values = pal)+
  ggtitle('PCA Iris dataset comparing 5 runs and 5 epsilons')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))


plot(g1)
plot(g2)
plot(g3)

pca<-prcomp(orig)
str(pca)
g4 <- ggplot(data.table(pca$x), aes(PC1, -PC2))+geom_point()+
  scale_color_manual(values = pal)+
  ggtitle('PCA Iris dataset comparing 5 runs and 5 epsilons')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))
plot(g4)
plot(g2)
ggbiplot(pca, var.axes = F)




dat <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_simulation_multiple.tsv')

drop_sign <- function(x){
  if(sign(x)==-1) {
    -1
  }
  else{
    1
  }
}
#multiplier<-sapply(dat$V1, function(x) drop_sign(x))
#dat$V1<-dat$V1 * multiplier
#dat$V2<-dat$V2 * multiplier

dat$index <- rep(1:150,500)

colnames(dat)<-c('V1', "V2", 'V3', "V4", "split", "re", "epsilon", "delta", "index")
g4 <- ggplot(dat[split==1], aes(V1, V2))+geom_point()+
  ggtitle('PCA')+
  ylab('PC2')+xlab('PC1')+theme(legend.title = element_blank(), plot.title = element_text(size = 28))+
  geom_mark_circle(aes(color=as.factor(index)), expand = unit( 2, "mm"))+facet_grid(~split)+theme(legend.position = 'none')
g4

data = datasets::iris
pp = prcomp(data[,1:4], center=T, scale.=T)
pca = fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_standalonePCA.tsv')
m = fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_2.tsv')
n = fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_3.tsv')

g4+
  geom_point(data = as.data.table(m), aes(x=V1,y=V2, col = '2')) 

ggplot()+  geom_point(data = as.data.table(m), aes(x=V1,y=V2, col = '3'))

scaled1 <-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_Scaled.tsv')

proj <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_2eigen.tsv')
eig1 <-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_2eigen.tsv')
eig2 <-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_simulation_eigen210.010.01.tsv')
eig3 <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_simulation_eigen200.010.01.tsv')

proj <- fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_2eigen.tsv')
p1 <-fread('/home/anne/Documents/featurecloud/featurecloud-test-data/federated_DP_PCA/IRIS_multi_no_noise_2.tsv')
p1
dat
scaled1<-as.matrix(scaled1)
eig1<-as.matrix(eig1)
eig2<-as.matrix(eig2) 
eig3<-as.matrix(eig3) 

d1 <- scaled1 %*% eig1
d2<- scaled1 %*% eig2
d3 <-scaled1 %*% eig3[,1:2]

ggplot()+geom_point(data=as.data.table(d1), aes(V1,V2, col='1'))+geom_point(data=as.data.table(d1), aes(V1,V2, col='2'))+
  geom_point(data = dat[V6==0], aes(V1,V2))+geom_point(data=as.data.table(d3), aes(V1,V2, col='3'))
