require(data.table)
require(ggplot2)

### This script visualises the execution times for single and multisite runs from the log files generated
###


setwd('/home/anne/Documents/featurecloud/results_final/split/')


single<-c('Power iteration', 'Single site PCA', 'Subspace iteration')
multi<-c('Unequal split proxy', 'Unequal split subspace iteration')
 
fne<-c()

times<-list()
for(d in list.dirs(recursive = F)){
  fl<-file.path(d, '0.5', 'time_log.tsv')
  if(file.exists(fl)){
    times[[d]]<-fread(fl)
  }else{
   fne<-c(fne, fl)
  }
  
}

time<-rbindlist(times)
single_site<-ggplot(time[V1 %in% single], aes(V1, V2))+
  geom_boxplot()+
  theme(axis.text.x = element_text())+
  ggtitle('Execution times for single site runs')
multi_site<-ggplot(time[V1 %in% multi & V2<1000], aes(V1, V2))+
  geom_boxplot()+
  theme(axis.text.x = element_text())+
  ggtitle('Execution times for multi site runs')+
  scale_x_discrete(labels=c('Unequal split proxy' = 'Single round', 'Unequal split subspace iteration' = 'Subspace iteration'))


ggsave(single_site, file='/home/anne/Documents/paper_fed_PCA/figures/execution_times_single_site.pdf')
ggsave(multi_site, file='/home/anne/Documents/paper_fed_PCA/figures/execution_times_multi_site.pdf')  

time$group<-as.numeric(as.factor(time$V1))
time$group<-ifelse(time$group %in% c(1,2,3,4), 'Global', 'Pooled')
all<-ggplot(time[V1 %in% c("Single site PCA",
                           "Subspace iteration",            
                           "Power iteration",
                           "Unequal split proxy" ,            
                           "Unequal split subspace iteration")], aes(V1, V2, fill = as.factor(group)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),)+
  ggtitle('Execution times')+
  scale_x_discrete(labels=c('Unequal split proxy' = 'Single round', 
                            'Unequal split subspace iteration' = 'Unequal split\nsubspace iteration'))+
  theme(axis.title.x = element_blank())+
  ylab('Wallclock time [s]')+scale_fill_discrete(name = "Model type")
all


nit <-list()
i = 1
for(set in list.dirs(recursive = F)){
  nr_iter<-  fread(file.path(set,'0.5/nr_iter_powit.tsv'))
  nr_iter <-data.table(V1 = nr_iter$V1, V2 = sum(nr_iter[, 2:21]), algo = 'Power\nIteration\n(Deflation)')
  nit[[i]]<-nr_iter
  i <- i+1
  nr_iter<-  fread(file.path(set,'0.5/nr_iter_qr.tsv'))
  nr_iter$algo <- 'Subspace\nIteration'
  nit[[i]]<-nr_iter
  i <- i+1
}
nit<-rbindlist(nit, fill=T)
colnames(nit)<-c('dataset', 'nr_iter', 'algo')


it <-ggplot(nit, aes(algo, nr_iter))+geom_boxplot()+
  labs(title = "Global models", y = '# Iterations', x = 'Algorithm')+
  theme(plot.title = element_text(size = 35, hjust = 0.5), 
        legend.key.size = unit(0.7, "cm"),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 20), 
        axis.text = element_text(size=20),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        legend.text=element_text(size=20))+
  labs(subtitle = 'Total # iterations required to retrieve all eigenvectors')+
  scale_fill_brewer('Model type', palette = 'Blues')
it
ggsave(it, file = '/home/anne/Documents/paper_fed_PCA/figures/number_iterations.pdf', height = 10, width = 9)
it
