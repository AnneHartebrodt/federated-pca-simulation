require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
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


#data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist10.encryption')
#data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist5.encryption')
#data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist.encryption')
data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist2.encryption')
summarised <- data %>% group_by(V1) %>% summarise(time = sum(V4))
summarised<- as.data.table((summarised))
summarised$time <- summarised$time/max(data$V2)
summarised$time<-summarised$time/2
summarised$time <- round(summarised$time, 2)
#summarised$time <- summarised$time/(sum(summarised$time))


selection<-c('matrix_encryption', 'matrix_addition', 'matrix_decryption', 'qr_encryption', 'qr_addition', 'qr_decryption')
byiteration<-summarised[V1 %in% selection]
colnames(byiteration)<-c('Step', 'Time [s]')
#fwrite(byiteration, file='/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist10.encryption.summary', sep='\t')
#fwrite(byiteration, file='/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist5.encryption.summary', sep='\t')
#fwrite(byiteration, file='/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist.encryption.summary', sep='\t')
fwrite(byiteration, file='/home/anne/Documents/featurecloud/pca/vertical-pca/results/encryption/mnist2.encryption.summary', sep='\t')