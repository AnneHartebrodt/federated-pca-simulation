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


encryption<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_encryption.encryption')
encryption2<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_encryption_encrypted.encryption')

encryption<-rbind(encryption, encryption2)
colnames(encryption)<-c("encrypted", 'sites', 'repeat', 'wallclock')

summary <- encryption %>% group_by(encrypted, sites) %>% summarise(mean_time = mean(wallclock))
summary <- as.data.table(summary)
#summary$mean_time<-summary$mean_time/summary$sites

summary$per_split <- summary$mean_time/ summary$sites 
summary$per_split <- round(summary$per_split, 3)
summary$mean_time<-round(summary$mean_time, 3)

# local rounds are averaged over the number of clients
# aggregation round scale with the number of clients and are therefore not averaged.
summary$time <- summary$per_split
summary[encrypted=='addition']$time <- summary[encrypted=='addition']$mean_time

fwrite(summary, file='/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/summary.txt', sep='&')

nr_sites = 5
barchart <- ggplot(summary[encrypted %in% c('qr_addition', 'qr_decryption', 'qr_encryption', 'qr_setup') & sites %in% c(nr_sites)], 
                   aes(x = factor(sites),y = time, fill=encrypted))+
  geom_bar(stat = 'identity')+
  scale_fill_manual('',values = palette_div[c(1,3,6,8)])+ ylab('Wallclock time [s]')+xlab('')+
  my_theme+coord_flip()+theme(legend.position = 'bottom', axis.ticks.y = element_blank(), axis.text.y = element_blank())
barchart
ggsave(barchart, file = '/home/anne/Documents/featurecloud/pca/vertical-pca/figures/qr/encryption_time.pdf', height=2)

barchart
s1 <- sum(summary[encrypted %in% c('addition', 'encryption', 'decryption', 'setup')]$mean_time)
s2 <- sum(summary[encrypted %in% c('total_encrypted')]$mean_time)
s3 <- sum(summary[encrypted %in% c('total_unencrypted')]$mean_time)

angles<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_angles_all_against_all.tsv')
angles<- angles %>% pivot_longer(-c(V1,V2, V3))
angles<-as.data.table(angles)
angles<-angles[!is.na(value)]
mean(angles$value)

angles_enc <- fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_angles_all_against_all_encrypted.tsv')
angles_enc<- angles_enc %>% pivot_longer(-c(V1,V2, V3))
angles_enc<-as.data.table(angles_enc)
angles_enc<-angles_enc[!is.na(value)]
mean(angles_enc$value)


enc_vs_unenc<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_angles_encrypted_vs_normal.tsv')
enc_vs_unenc<- enc_vs_unenc %>% pivot_longer(-c(V1,V2, V3))
enc_vs_unenc<-as.data.table(enc_vs_unenc)
enc_vs_unenc<- enc_vs_unenc[!is.na(value)]
sum(enc_vs_unenc$value)


enc_vs_unenc<-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_angles_scipy_vs_federated.tsv')
enc_vs_unenc<- enc_vs_unenc %>% pivot_longer(-c(V1,V2, V3))
enc_vs_unenc<-as.data.table(enc_vs_unenc)
enc_vs_unenc<- enc_vs_unenc[!is.na(value)]
sum(enc_vs_unenc$value)/ncol(enc_vs_unenc)

transmission <-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_encryption.transmission')
transmission.encrpted <-fread('/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption/benchmark_encryption_qr_encryption_encrypted.transmission')

