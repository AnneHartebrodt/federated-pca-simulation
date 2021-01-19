require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(plyr)
require(dplyr)
require(R.utils)
library(ggthemes)   
require(scales)
require(R.utils)
library(ggforce)
require(stringr)


palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.key.size = unit(1.5, 'lines'))



read_directory<-function(current.directory, file.ending, h,m,s,q,g){
  angle_list<-list()
  counter<-1
  myfiles <-  list.files(file.path(current.directory), recursive = F)
  for (f in which(str_detect(myfiles, file.ending))) {
    filename <- myfiles[f]
    #print(filename)
    #print(file.path(current.directory, filename))
    if(!file.exists(file.path(current.directory, filename))){
      next
    }
    if(file.size(file.path(current.directory, filename))==0){
      next
    }
    angles <-fread( file.path(current.directory, filename), sep = '\t',header = F)
    angles<- angles[,1:(ncol(angles)-1)]
    
    if(ncol(angles)>3){
    # transform data into correct format
    colnames(angles) <-c('iterations', sapply(1:(ncol(angles) - 1), function(x) paste0('', x)))
    
    angles$filename <- filename
    angles$partitioning<-h
    angles$split <- s
    angles$qr<-q
    angles$evupdate<-g
    angles$method<-m
    angles <-as.data.table(pivot_longer(angles,-c(iterations,partitioning, split, qr, filename, evupdate, method)))
    }else{
    colnames(angles)<-c( 'name', 'iterations', 'value')
    angles<-angles[,c('iterations', 'name', 'value')]
    angles$filename <- filename
    angles$partitioning<-h
    angles$split <- s
    angles$qr<-q
    angles$evupdate<-g
    angles$method<-m
  }
    
    angle_list[[counter]]<-angles
    counter<-counter+1
  }
  angles<-rbindlist(angle_list)
  angles<-as.data.frame(angles)
  if(nrow(angles)!=0){
  angles<-angles[, c("iterations","partitioning","split","qr","evupdate","name","method", "value")]
  return(angles)
  }
  else{
    return(NA)
  }
}

read_data<-function(basedir, file.ending){
counter = 1
angle_list<-list()
for(h in c('horizontal', 'vertical')){
  for(m in  c('matrix', 'vector')){
  for(s in c(10, 2, 3, 5)){
    for(g in c('power', 'gradient')){
      for(q in c('central_qr', 'federated_qr')){
        current.directory <- file.path(basedir, h,m,s,g, q)
        #print(current.directory)
        # pass if directory does not exist
        if(!file.exists(current.directory)){
          next
        }
        angles<-read_directory(current.directory , file.ending,h,m,s,q,g)
        angle_list[[counter]]<-angles
        counter <- counter+1
        }
      }
    }
  }
  for (m in c("power_iteration")){
    for(s in c(10, 2, 3, 5)){
    current.directory <- file.path(basedir,h,m,s)
    # pass if directory does not exist
    if(!file.exists(current.directory)){
      next
    }
    angles<-read_directory(current.directory , file.ending, h,m,s,NA, NA)
    angle_list[[counter]]<-angles
    counter <- counter+1
    }
  }
}


angles<-rbindlist(angle_list)
angles<-as.data.table(angles)
angles$name<-as.numeric(angles$name)
return(angles)
}

basedir<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist_test'

angles<-read_data(basedir, file.ending = 'angles.u')
angles$qr<-revalue(angles$qr, c("central_qr"="CENT", "federated_qr"="FED"))
angles$evupdate<-revalue(angles$evupdate, c("gradient"="GRAD", "power"="POW"))
angles$method<-revalue(angles$method, c("matrix"="SIM", "vector"="SEQ"))
angle.summary<-angles %>% group_by(iterations, partitioning, split, qr, evupdate, name, method) %>% summarise(mean_angle=mean(value))
angle.summary<-as.data.table(angle.summary)
angle.summary$name<-as.numeric(angle.summary$name)
angle.summary[,Experiment:=paste(method, qr, evupdate, sep='-')]

for( s in c(2,3,5,10)){
s2<-ggplot(angle.summary[split==s & partitioning=='horizontal'], aes(iterations, mean_angle, col=Experiment))+
  geom_line()+
  my_theme+
  scale_color_manual(values = palette_div)+
  ggtitle('Horizontal partitioning - Right Singular Vectors')+
  xlim(0,50)+
  xlab('#Iterations')+
  ylab('Mean angle [degree]')+
  facet_wrap(~name+method)
ggsave(s2, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/manuscript/mnist_hor_',s ,'left.pdf'), width = 30, units = 'cm', height = 20)
}

for( s in c(2,3,5,10)){
ggplot(angle.summary[split==s & partitioning=='vertical'], aes(iterations, mean_angle, col=Experiment))+
  geom_line()+
  my_theme+
  scale_color_manual(values = palette_div)+
  ggtitle('Vertical - Right Singular Vectors')+
  xlim(0,50)+
  xlab('#Iterations')+
  ylab('Mean angle [degree]')+
  facet_wrap(~name+method)
  ggsave(s2, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/manuscript/mnist_ver_',s ,'right.pdf'), width = 30, units = 'cm', height = 20)
}


angles.v<-read_data(basedir, file.ending = 'angles.v')
angles.v$qr<-revalue(angles.v$qr, c("central_qr"="CENT", "federated_qr"="FED"))
angles.v$qr[is.na(angles.v$qr)]<-''
angles.v$evupdate<-revalue(angles.v$evupdate, c("gradient"="GRAD", "power"="POW"))
angles.v$evupdate[is.na(angles.v$evupdate)]<-''
angles.v$method<-revalue(angles.v$method, c("matrix"="SIM", "vector"="SEQ", "power_iteration"="POW-HOR"))
angle.summary.v<-angles.v %>% group_by(iterations, partitioning, split, qr, evupdate, name, method) %>% summarise(mean_angle=mean(value))
angle.summary.v<-as.data.table(angle.summary.v)
angle.summary.v$name<-as.factor(as.numeric(angle.summary.v$name))
angle.summary.v[,Experiment:=paste(method, qr, evupdate, sep='-')]
angle.summary.v$Experiment[angle.summary.v$Experiment=='POW-HOR--']<-'POW-HOR'

for( s in c(2,3,5,10)){
s2<-ggplot(angle.summary.v[split==s & partitioning=='horizontal'], aes(iterations, mean_angle, col=Experiment))+
  geom_line()+
  my_theme+
  scale_color_manual(values = palette_div)+
  ggtitle('Horizontal partitioning - Left Singular Vector')+
  xlim(0,50)+
  xlab('#Iterations')+
  ylab('Mean angle [degree]')+
  facet_wrap(~name+method)
ggsave(s2, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/manuscript/mnist_hor_',s ,'right.pdf'), width = 30, units = 'cm', height = 20)
}

  


for( s in c(2,3,5,10)){
ggplot(angle.summary.v[split==s & partitioning=='vertical'], aes(iterations, mean_angle, col=Experiment))+
  geom_line()+
  my_theme+
  scale_color_manual(values = palette_div)+
  ggtitle('Vertical partitioning - Left Singular Vector')+
  xlim(0,50)+
  xlab('#Iterations')+
  ylab('Mean angle [degree]')+
  facet_wrap(~name+method)
ggsave(s2, file=paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/manuscript/mnist_ver_',s ,'left.pdf'), width = 30, units = 'cm', height = 20)
}