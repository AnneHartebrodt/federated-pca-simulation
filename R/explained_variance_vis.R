require(data.table)
require(ggplot2)

setwd('/home/anne/Documents/featurecloud/')
variance<-fread('results/pca_plots/all/var_explained_aor.txt')
colnames(variance)<-as.character(seq(0.1, 1, by=0.1))
studies<-read.table('results/pca_plots/all/study_names.tsv')
variance$study.id<-studies$V1
variance<-melt(variance, value.name = 'nr.vars', id.vars = 'study.id')


study.sizes<-fread('/home/anne/Documents/featurecloud/results/usability_study/TCGA_study_sizes.tsv')
colnames(study.sizes)<-c('study.id', 'nr.samples', 'database')
study.sizes[study.id=='BEATAML1.0-COHORT', 1]<-'BEATAML1'

variance<-merge(variance, study.sizes, by = 'study.id')
variance[, perc:= nr.vars/nr.samples]

ggplot(variance[variable %in% c('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8')], aes(x = variable,y = perc))+
  geom_boxplot()+
  ylab('#PCs/#PCs to explain 100% of variance')+
  xlab('Explained variance')

ggplot(variance[variable %in% c('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8')], aes(x = variable,y =nr.vars))+
  geom_boxplot()+
  ylab('#PCs/#PCs to explain 100% of variance')+
  xlab('Explained variance')

ggplot(variance, aes(x = variable,y = perc))+
  geom_boxplot()+
  ylab('#PCs/#PCs to explain 100% of variance')+
  xlab('Explained variance')

ggplot(variance, aes(x = variable,y =nr.vars))+
  geom_boxplot()+
  ylab('#PCs/#PCs to explain 100% of variance')+
  xlab('Explained variance')
