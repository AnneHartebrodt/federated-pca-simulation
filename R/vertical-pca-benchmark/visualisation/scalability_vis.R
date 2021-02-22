require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)

dat<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/results_for_david//wide.tikz.iterations.transmission_cost.tsv')

da<- dat %>% pivot_longer(-counter)
da<-as.data.table(da)


ggplot(da, aes(name, value))+geom_boxplot()+theme(axis.text.x = element_text(angle=90))
