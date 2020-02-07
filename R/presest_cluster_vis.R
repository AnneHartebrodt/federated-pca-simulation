require(data.table)
require(ggplot2)

dir<-'/home/anne/Documents/featurecloud/results/split_study_cluster'

data<-NULL
annot<-data.table(var=c(0.5, 0.5, 0.2, 0.2, 0.8,0.8), clu =c(2,3,2,3,2,3))
for ( d in list.files(dir)){
  if(file.exists(file.path(dir, d, 'angles_cluster_splits.tsv'))){

      dda<-fread(file.path(dir, d, 'angles_cluster_splits.tsv'))
      dda <- cbind(annot, dda)
      data[[d]]<-dda
  }
}
data<-rbindlist(data)
data<-melt(data, id.vars = c('var', 'clu', 'V1'))
ggplot(data[var==0.5], aes(x = variable, y=value, fill = as.factor(clu)))+geom_boxplot()#+facet_wrap(~var)
