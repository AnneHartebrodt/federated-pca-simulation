require(data.table)
require(ggplot2)
require(tidyr)
require(stringr)
source('/home/anne/Documents/featurecloud/pca/federated_dp_pca/R/horizontal-pca-benchmark/library.R')

type = "Corpus_uteri"

data<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.tsv'))
names<-fread(paste0('/home/anne/Documents/featurecloud/data/tcga/cancer_type/',type,'/sites/data_all/',type,'.names.tsv'),
             header=F, sep='\t')

names[,c('type', 'tss', 'id') := tstrsplit(V1, '-')]


pca <- prcomp(data, scale. = F, center = F)

x <- as.data.table(pca$x)[, 1:10]
x$id<-names$id
x$tss<-names$tss
x$name<-names$V1


ggplot(x, aes(PC1, PC2, col=as.factor(tss), shape=tss))+
  geom_point(aes(shape=tss))+
  my_theme+
  scale_color_manual(name='TSS',values = c(palette_div, palette_div))+
  scale_shape_manual(name = 'TSS',values = c(1:20))+
  guides(color=guide_legend(ncol=2))+
  geom_label_repel(data = subset(x, name %in% ol),
                   aes(label=name), size=3)


long<-list()
for (p in c(colnames(x)[1:10])){
  y <- x %>% pivot_longer(-c(name, p), names_to = 'PC_B')
  y$PC_A <- p
  colnames(y)<-c('value.A', 'name', 'PC_B', 'value.B','PC_A')
  long[[p]]<-y
}
long<-rbindlist(long)


ggplot(long, aes(value.A, value.B, col=name, shape=name))+
  geom_point(aes(shape=name))+
  facet_grid(c('PC_B', 'PC_A'), scales = 'free')+
  scale_color_manual(name='TSS',values = c(palette_div, palette_div))+
  scale_shape_manual(name = 'TSS',values = c(1:20))+
  guides(color=guide_legend(ncol=2))+
  ylab("")+xlab("")



ggplot(x, aes(PC1, PC3, col=as.factor(name), shape=name))+
  geom_point(aes(shape=name))+
  my_theme+
  scale_color_manual(name='TSS',values = c(palette_div, palette_div))+
  scale_shape_manual(name = 'TSS',values = c(1:20))+
  guides(color=guide_legend(ncol=2))


ggplot(as.data.table(pca$x), aes(PC2, PC3, col=as.factor(names$V1), shape=names$V1))+
  geom_point(aes(shape=names$V1))+
  my_theme+
  scale_color_manual(name='TSS',values = c(palette_div, palette_div))+
  scale_shape_manual(name = 'TSS',values = c(1:20))+
  guides(color=guide_legend(ncol=2))
