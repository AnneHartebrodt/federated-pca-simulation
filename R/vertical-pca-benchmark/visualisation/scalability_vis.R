require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(viridis)
library(ggpubr)
require(optparse)

option_list = list(
  make_option(c("-f", "--infile"), action="store", default=NA, type='character',
              help="infile"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="outfile")
  
  
)
opt = parse_args(OptionParser(option_list=option_list))

infile<-opt$infile
infile <- '/home/anne/Documents/featurecloud/pca/vertical-pca/results-new-tests/scalability/summary.transmission_cost.tsv'
outfile<-opt$outfile
outfile<- '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/transmission_cost.pdf'
column<-opt$column
print(infile)


data<-fread(infile)

data$cost<- data$nr_float*4/(1000000)
#data$facet_title <- ordered(data$facet_title, levels=levels(data$rank))
data[, name_qr := paste0(matrix, eigenvector_update, qr_method, orthonormalisation_skip)]
data$rank <-as.factor(data$rank)
data$name_qr<-as.factor(data$name_qr)
data$name_qr<- recode(data$name_qr, matrixpowerfederated_qr1='FED-GS', matrixpowerfederated_qr100='NO-GS', vectorgradientcentral_qr1='GUO')
data <- as.data.table(data)



cost<- ggplot(data, aes(as.factor(data_fraction), cost, fill=name_qr))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.line=element_line(),
      legend.position = 'right',
      strip.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size=6), 
      legend.title = element_blank(),
      title = element_text(size=8),
      legend.text = element_text(size = 6),
      ,
      text=element_text(family = 'Arial'))+
  ylab('transmitted data volume in MB')+
  xlab('percentage of samples')+
  ggtitle('Estimated size of transmitted data')+
  scale_fill_manual(values = viridis(4)[1:3])+ guides(color=guide_legend(keyheight = 0.6, title = element_text('Configuration', size = 8)))
cost

ggsave(cost, filename = outfile , height = 6, units = 'cm', width = 10, device = cairo_pdf)

iterations <-ggplot(data, aes(as.factor(data_fraction), max_it, fill=name_qr))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.line=element_line(),
        legend.position = 'right',
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size=7), 
        legend.title = element_blank())+
  ylab('iterations')+
  xlab('percentage of samples')+
  ggtitle('Total number of iterations')+
  scale_fill_manual(values = viridis(4)[1:3])

