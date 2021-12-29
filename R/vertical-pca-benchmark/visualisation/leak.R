require(data.table)
require(ggplot2)
require(tidyverse)
data<-fread('/home/anne/Documents/featurecloud/pca/approximative-vertical/results/leeka.tsv')
colnames(data)<-c('counter','k', 'features', 'nsamples','restricted', 'maxit', 'type', 'error', 'time')

dat


extra_theme <-theme(axis.text.y = element_text( size=8),
                    axis.title.y = element_text(size=9))
                    
gp<-ggplot(data[k %in% c(2,5) & type %in% c('reconstructed', 'baseline')], aes(as.factor(features), error, color=as.factor(type), fill=as.factor(nsamples), group=interaction(features, nsamples, type)))+
  geom_boxplot()+
  theme_classic()+
  xlab('#Features')+
  facet_wrap('k')+
  ylab('Correlation of covariance matrices')+extra_theme+
  scale_fill_viridis_d('#Samples',option='D', end=0.7)+
  scale_color_viridis_d('Reconstruction\ntechnique',option='A', end=0.7)
gp
ggsave(gp,file='/home/anne/Documents/manuscripts/vertical-pattern-recognition/figures/leak.pdf', height = 8, width=20, unit='cm')

