require(data.table)
require(ggplot2)

dir<-"~/Documents/featurecloud/results/split_study_cluster/CPTAC-2/"

meta<-fread(paste0(dir, 'meta_splits.tsv'), fill = T)
combo<-unique(meta[,2:(ncol(meta)-1)])

group<-factor(rep(c('1.0', 
                '0.1/0.9', '0.3/0.7', '0.5/0.5',
                '0.2/0.2/0.2/0.2/0.2', '0.1/0.1/0.2/0.2/0.4', '0.1/0.1/0.1/0.1/0.6', '0.2375/0.2375/0.2375/0.2375/0.05'), each=10))
angs = list()
i  = 1
for( d in c('1.0', '1.5', '2.0', '5.0')){
  ang1<-fread(paste0(dir, 'angles_unequal_splits',d,'.tsv'), fill = T)
  ang1$group<-group
  ang1$split<-d
  angs[[i]]<-ang1
  i = i+1
}


ang<-rbindlist(angs)
ang$split<-as.factor(ang$split)
ggplot(ang, aes(group, V2))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90))

