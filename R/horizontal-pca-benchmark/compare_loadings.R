require(data.table)
require(ggplot2)
require(colorspace)
require(RColorBrewer)
# make custom palettes
palette_div<-c('#062865', '#2f497d', '#4c6d96', '#6793af', '#81bac8', '#ffc4b3', '#f68888', '#d75161', '#ab203f', '#720022')
palette_seq<-c('#062865', '#203a72', '#324d80', '#43618d', '#52759b', '#618aa9', '#70a0b7', '#7eb6c5', '#8dcdd4', '#9ce4e2')

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.key.size = unit(1.5, 'lines'))

top<-fread('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/single-cell/overlap_top.tsv', header=T)
low<-fread('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/single-cell/overlap_bottom.tsv', header=T)
low.paga <-fread('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/single-cell/paga_overlap_bottom.tsv', header=T)
top.paga<-fread('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/single-cell/paga_overlap_top.tsv', header=T)

require(dplyr)
require(tidyr)
df<-data.table(count = 1:nrow(top), PBMC.top=top$`0`, PBMC.low=low$`0`, PAUL.top = top.paga$`0`, PAUL.low = low.paga$`0`)
df<- as.data.table(df %>% pivot_longer(-count))

colnames(df)<- c('count', 'name', 'Genes')
g1<-ggplot(df[count<11], aes( name, as.factor(count),size=Genes, color=Genes, label=Genes))+
  geom_point()+
  ylab("Eigenvector rank")+
  scale_colour_gradient(low= '#062865', high ='#81bac8',
                        guide = guide_legend(
                          direction = "vertical",
                          title.position = "bottom",
                          label.position = "right",
                          label.vjust = 0.5,
                          label.theme = element_text(size=8, vjust = -1),
                          keywidth = unit(0.5, 'cm'), keyheight =unit(0.5, 'cm'),
                          title.theme = element_text(size = 10, hjust = 0.5))) +
  scale_size_continuous()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5), 
        axis.title.x = element_blank())+
  scale_y_discrete(limits=rev)+
  scale_x_discrete(limits=rev)

g1

ggsave(g1, file='/home/anne/Documents/manuscripts/horizontal-pca/figures/genes_overlap.pdf', height = 7.6, width = 6, units = 'cm')
