require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(tidyverse)



data1<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix/summary.angles.u.tsv')
data1$index<-1:nrow(data1)
colnames(data1)<- c('iterations', paste0('EV', 1:10), 'sites', 'algorithm', 'maxit', 'filename', 'index')

# get last iteration and index
max_rows <- data1 %>% group_by( sites, algorithm, maxit, filename)  %>% slice_max(iterations)
iterations<-max_rows %>%ungroup() %>% select(iterations, sites, algorithm)

globalmax<- data1[,max(iterations)]
repeats<- abs(max_rows$iterations-globalmax)
index<-unlist(sapply(max_rows$iterations[max_rows$iterations!=globalmax], function(x) seq(from=x+1, globalmax)))
max_rows <- data1[rep(x = max_rows$index, repeats)]

max_rows$iterations<-index
max_rows<- as.data.table(max_rows)
data1 <- rbind(data1, max_rows)

summary1 <- data1  %>% group_by(iterations, sites,algorithm) %>% summarise(across(c(paste0('EV', 1:10)),mean)) %>% pivot_longer(-c(iterations, sites,algorithm))
summary1<-as.data.table(summary1)
summary1$algorithm<-as.factor(summary1$algorithm)
colnames(summary1)<-c('iterations', 'sites', 'algorithm', 'eigenvector', 'value')
ggplot(summary1, aes(iterations, value,color=algorithm))+geom_line(aes(linetype=algorithm))+facet_grid(c( 'eigenvector','sites' ))



data<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.angles.u.tsv')
colnames(data)<-  c('eigenvector', 'iterations', 'value', 'sites','algorithm', 'dummy', 'filename')
max_it<- data %>% group_by( sites, algorithm, filename)  %>% slice_max(iterations) %>% ungroup() %>%  select(iterations, sites, algorithm)
data<- data  %>%group_by(eigenvector, filename) %>% mutate(counter = row_number(iterations))
offset<-data %>% group_by(eigenvector, sites, algorithm, filename)%>% summarise(offset=max(iterations)) %>% summarise(mean_offset = floor(mean(offset)))
offset<- as.data.table(offset)
offset$mean_offset<-offset$mean_offset+1
offset<-rbind(offset, list(eigenvector=1,sites=5, algorithm='guo', mean_offset = 0))
offset<-rbind(offset, list(eigenvector=1,sites=10, algorithm='guo', mean_offset = 0))
data <- data %>% group_by(counter, eigenvector, sites, algorithm)%>% summarise(avg_angle = mean(value))
data<-data %>% left_join(offset)


data<-as.data.table(data)

do<-list()
counter <-1
offset<-offset[order(V1)]
for (i in 3:20){
  print(i)
  do[[counter]]<-data.table(c(1:as.numeric(offset[i,4])), offset[i,1], offset[i,2], 'guo', 90, '1')
  counter<-counter+1
}


do<-rbindlist(do)

data$counter<-data$counter +data$mean_offset
colnames(do)<-colnames(data)

data<-rbind(data, do)
colnames(data)<-c('iterations', 'eigenvector', 'sites', 'algorithm', 'value', 'maxit')
data$eigenvector<-paste0('EV', data$eigenvector)
data <-data %>% select(iterations, eigenvector, sites, algorithm, value)
summary1<-rbind(summary1, data)
summary1$algorithm<-recode(summary1$algorithm, 'approximative-init'='AI-FULL', 
                           'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 
                           'randomized-approx-projected'='AI-RAND',
                           'guo'='GUO')
summary1 <- summary1[algorithm %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO')]
summary1$eigenvector<-ordered(summary1$eigenvector,  levels = paste0('EV', 1:10))
ggplot(summary1, aes(iterations, value,color=algorithm))+geom_line( size=1)+
  facet_grid(c( 'eigenvector','sites' ))+xlim(0,250)+
  xlab('Angles w.r.t reference')+theme_classic()+ylab('Iterations')+scale_color_viridis_d(option = 'A', end=0.7)



time <- fread('~/Documents/featurecloud/pca/approximative-vertical/results/time.log')
time$V1<-recode(time$V1, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
time<- time[V1 %in%  c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND', 'GUO')]
ggplot(time, aes(V1, V4))+geom_boxplot()+ylab('Runtime[s]')+xlab('')+theme_classic()+scale_color_viridis_d(option = 'A', end=0.7)


eo<- fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix/summary.eo.tsv')
eo$V4<-recode(eo$V4, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND')
eo<- eo[V4 %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND') ]
eo <- eo %>% group_by(V3,V4, V5, V6) %>% summarise(sum=sum(V2))
ggplot(eo, aes(V4, sum))+geom_boxplot()+ylab('Runtime (Matrix operations) [s]')+xlab('')+theme_classic()

iterations<-rbind(as.data.table(iterations), as.data.table(max_it))
iterations$algorithm<-recode(iterations$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
iterations <- iterations[algorithm %in% c('AI-FULL', 'RI-FULL', 'AI-RAND', 'RI-RAND') ]
ggplot(iterations, aes(algorithm, iterations, fill=as.factor(sites)))+geom_boxplot()+xlab('')+theme_classic()+ylab('Iterations until convergence')


mev<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/matrix//summary.mev.u.tsv')
colnames(mev)<-c('iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
mev<-as.data.table(mev)
mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
mev<- mev[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND', 'AI-FULL') ]
ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
  xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)


mev<-fread('~/Documents/featurecloud/pca/approximative-vertical/results/vector/summary.mev.u.tsv')
colnames(mev)<-c('vector','iterations', 'mev', 'sites', 'algorithm', 'maxit', 'filename')
mev <- mev %>% group_by(iterations, sites, algorithm) %>% summarise(mean_mev = mean(mev))
mev<-as.data.table(mev)
mev$algorithm<-recode(mev$algorithm, 'approximative-init'='AI-FULL', 'random-init'='RI-FULL', 'randomized-no-approx'='RI-RAND', 'randomized-approx-projected'='AI-RAND', 'guo'='GUO')
mev<- mev[algorithm %in% c( 'RI-FULL', 'AI-RAND', 'RI-RAND') ]
ggplot(mev, aes(iterations,mean_mev, color=algorithm))+geom_line()+facet_wrap(~sites)+
  xlab('MEV')+theme_classic()+ylab('iterations')+xlim(0,100)+scale_color_viridis_d(option = 'A', end=0.7)
