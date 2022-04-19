require(data.table)
require(ggplot2)
require(dplyr)

byte2giga <- 100000000

data<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/matrix_encryption/encryption_time.tsv', header=F)
summary <- data %>% group_by(V2, V4) %>% summarise(mean_time = mean(V5), mean_size_ne =mean(V6), mean_size_enc=mean(V7)) 
summary <- as.data.table(summary)
summary$mean_size_enc<- summary$mean_size_enc/byte2giga
summary$mean_size_ne <- summary$mean_size_ne/byte2giga

ggplot(data[], aes(as.factor(V2),V5))+geom_boxplot()+scale_y_log10()+facet_wrap(~V4, scales = 'free')+
  xlab('Elements in matrix')+ylab('Encryption time [s]')
    