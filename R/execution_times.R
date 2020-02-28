require(data.table)
require(ggplot2)

setwd('/home/anne/Documents/featurecloud/results/sandbox/splits/')


single<-c('Power iteration', 'Single site PCA', 'Subspace iteration')
multi<-c('Unequal split proxy', 'Unequal split subspace iteration')

times<-list()
for(d in list.dirs(recursive = F)){
  times[[d]]<-fread(file.path(d, '0.5', 'time_log.tsv'))
}

time<-rbindlist(times)
single_site<-ggplot(time[V1 %in% single], aes(V1, V2))+
  geom_boxplot()+
  theme(axis.text.x = element_text())+
  ggtitle('Execution times for single site runs')
multi_site<-ggplot(time[V1 %in% multi], aes(V1, V2))+
  geom_boxplot()+
  theme(axis.text.x = element_text())+
  ggtitle('Execution times for multi site runs')+
  scale_x_discrete(labels=c('Unequal split proxy' = 'Single round', 'Unequal split subspace iteration' = 'Subspace iteration'))


ggsave(single_site, file='/home/anne/Documents/paper_fed_PCA/figures/execution_times_single_site.pdf')
ggsave(multi_site, file='/home/anne/Documents/paper_fed_PCA/figures/execution_times_multi_site.pdf')  
