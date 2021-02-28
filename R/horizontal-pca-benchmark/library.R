require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require("ggrepel")

# make custom palettes
palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')

# create a minimal theme
my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.key.size = unit(1.5, 'lines'))

read_cancer_types<-function(outdir){
  dirs<-list.files(outdir, no..=F)
  data_list<-list()
  for (d in dirs){
    data<-fread(file.path(outdir,d, paste0(d, '_angles.tsv')), sep='\t', fill=TRUE)
    
    colnames(data)<- c('algorithm', 1:(ncol(data)-1))
    data$dataset<-d
    data<- data %>% pivot_longer(-c(algorithm,dataset))
    data<-as.data.table(data)
    data_list[[d]]<-data
  }
  data <- rbindlist(data_list)
  data<-data[!is.na(value)]
  data$name<-as.factor(as.numeric(data$name))
  return(data)
}