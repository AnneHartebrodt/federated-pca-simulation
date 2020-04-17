inputfile = '/home/anne/Documents/featurecloud/results/sandbox2/TCGA-KIRC_TCGA-LIHC.sub.tsv'
outfile = '/home/anne/Documents/featurecloud/results/sandbox2'
clusterfile = '/home/anne/Documents/featurecloud/results/sandbox2/TCGA-KIRC_TCGA-LIHC_clusters.tsv'

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  t<-dot.prod / (norm.x * norm.y)
  #print(t)
  #Make the angle save to numerical inaccuracies
  if(t<=(-1)){
    #print('t<1')
    t<--1
  }
  if(t>1){
    #print('t>1')
    t<-1
  }
  theta <- acos(t)
  #print(theta)
  #print(as.numeric(theta))
  a<-as.numeric(theta)*180/pi
  #print(a)
  if(is.nan(a)){
    #print(a)
    return(1)
  }
  #print(a)
  a
}

require(ggplot2)
require(data.table)
require(gridExtra)

data <-fread(inputfile)
data <-data+1
data<-log2(data)
data<-scale(data, scale = FALSE, center = TRUE)

datrand <-data[sample(1:nrow(data), size = 400, replace = F),]

custom_cov <-function(data){
  cc <- (1/(nrow(data)-1))*(t(data)%*% data)
  return (cc)
}

cluster <- fread(clusterfile)

p = prcomp(data,  center = FALSE, scale = F)
pcaGlob <-prcomp(cov(data), center = FALSE, scale = F)
angle(p$rotation[,1],pcaGlob$rotation[,1])
Uall = pcaGlob$rotation[,1:30] %*% diag(x = p1$sdev[1:30], 30,30)


prand = prcomp(datrand,  center = FALSE, scale = F)
pcaGlobR <-prcomp(cov(datrand), center = FALSE, scale = F)
angle(p$rotation[,1],pcaGlob$rotation[,1])
UallR = pcaGlobR$rotation[,1:30] %*% diag(x = p1$sdev[1:30], 30,30)

data1<-data[which(cluster$V2 ==1),]
data2 <- data[cluster$V2 == 2,]
l1 = nrow(data1)/nrow(data)
l2 = nrow(data2)/nrow(data)


p1 = prcomp(data1,  center = FALSE, scale = F)
pcaGlob1 <-prcomp(custom_cov(data1), center = FALSE, scale = F)
angle(p1$rotation[,1],pcaGlob1$rotation[,1])

p2 = prcomp(data2,  center = FALSE, scale = F)
pcaGlob2 <-prcomp(custom_cov(data2), center = FALSE, scale = F)
angle(p2$rotation[,1],pcaGlob2$rotation[,1])

U = p1$rotation[,1:30] %*% diag(x = p1$sdev[1:30], 30,30)
prox1 <-  U %*% t(U)
U2  = p2$rotation[,1:30] %*% diag(x = p2$sdev[1:30], 30,30)
prox2 <- U2  %*% t(U2)

prox <- (l1*prox1+l2*prox2)
pcaPool <- prcomp(prox, center = F, scale = F)
Upool = pcaPool$rotation[,1:30] %*% diag(x = p1$sdev[1:30], 30,30)
t = angle(pcaPool$rotation[,1],pcaGlob$rotation[,1])


global_cov = custom_cov(data)
di = abs(prox - global_cov)
sum(di)
sum(abs(global_cov))
cd = melt(di)

cd[abs(value)>1]

cd$log = log2(abs(cd$value))
cd = as.data.table(cd)
t = cd[abs(value)>0.5]
cd[abs(value)<0.1]$value<-NA
ggplot(cd, aes(log))+geom_histogram(bins = 100)

ggplot(cd, aes(Var1, Var2, fill=value)) + geom_tile()+theme(axis.text = element_blank())

da = data.table(PC1 = U[,1], PC2 = U2[,1], f = (abs(U[,1])-abs(U2[,1])))
ggplot(da, aes(PC1, PC2, color = f))+geom_point()+geom_abline(slope = -1)+
  geom_abline(slope = -1, intercept = 0.8)+geom_abline(slope = -1,intercept = -0.8)

uo <-sort(U[,1], decreasing = T)
uo2 <-sort(U[,1], decreasing = T)
d<-intersect(names(uo[1:100]), names(uo2[1:100]))
uu2 = U2[which(rownames(U2) %in% d),]
uu = U[which(rownames(U) %in% d),]
uu<-melt(uu)
u2<-melt(uu2)

uuall = Uall[which(rownames(Uall) %in% d),]
uuall <- melt(uuall)

upool = Upool[which(rownames(Upool) %in% d),]
upool <- melt(upool)

uallR = UallR[which(rownames(UallR) %in% d),]
uallR <- melt(uallR)

uallR$diff = uallR$value-upool$value
uuall$diff = uuall$value -upool$value
r1<-ggplot(uu, aes(Var1, Var2, fill =value))+geom_tile()+theme(axis.text = element_blank())+ggtitle('prox1')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r2<-ggplot(u2, aes(Var1, Var2, fill =value))+geom_tile()+theme(axis.text = element_blank())+ggtitle('prox2')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r3<-ggplot(uuall, aes(Var1, Var2, fill =value))+geom_tile()+theme(axis.text = element_blank())+ggtitle('all')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r4<-ggplot(upool, aes(Var1, Var2, fill =value))+geom_tile()+theme(axis.text = element_blank())+ggtitle('pool')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r5<-ggplot(uallR, aes(Var1, Var2, fill =value))+geom_tile()+theme(axis.text = element_blank())+ggtitle('rnad')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r6<-ggplot(uallR, aes(Var1, Var2, fill =diff))+geom_tile()+theme(axis.text = element_blank())+ggtitle('rand-pool')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
r7<-ggplot(uuall, aes(Var1, Var2, fill =diff))+geom_tile()+theme(axis.text = element_blank())+ggtitle('all-pool')+scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") 
U = U[order(U[,1]),]
ggplot(data.table(y = abs(-U[,1]), x=1:length(U[,1]), t = '1'), aes(x,y, color = t))+geom_point()+geom_point(data = data.table(y = abs(-sort(U2[,1])), x=1:length(U[,1]),t = '2'), aes(x,y, color=t))+
  geom_point(data = data.table(y = abs(-uallRo[,1]), x=1:length(uallRo[,1]), t = '3'), aes(x,y, color = t))+
  geom_point(data = data.table(y = abs(sort(Uall[,1])), x=1:length(Uall[,1]), t = '4'), aes(x,y, color = t))

uallRo = UallR[order(UallR[,1]),]
ggplot(data.table(y = abs(uallRo[,1]), x=1:length(uallRo[,1])), aes(x,y))+geom_point()

grid.arrange(r1,r2,r3,r4, r5, r6, r7)
