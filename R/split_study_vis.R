require(data.table)
require(ggplot2)

dir<-"~/Documents/featurecloud/results/gexp_stats/target/BEATAML1/0.52.075396454/"

meta<-fread(paste0(dir, 'meta_splits.tsv'), fill = T)
combo<-unique(meta[,2:(ncol(meta)-1)])

group<-factor(rep(c( 
                '0.1/0.9', '0.3/0.7', '0.5/0.5',
                '0.2/0.2/0.2/0.2/0.2', '0.1/0.1/0.2/0.2/0.4', '0.1/0.1/0.1/0.1/0.6', '0.2375/0.2375/0.2375/0.2375/0.05'), each=10))
angs = list()
i  = 1
for( d in c('0.75', '10.0', '5.0')){
  ang1<-fread(paste0(dir, 'angles_unequal_splits_weighted_',d,'.tsv'), fill = T)
  ang1$group<-group[1:nrow(ang1)]
  ang1$split<-d
  angs[[i]]<-ang1
  i = i+1
}


ang<-rbindlist(angs, fill = T)
ang$split<-as.factor(ang$split)

ggplot(ang, aes(group, V2))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90))


ggsave(ggplot(ang, aes(group, V2))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig1.pdf')
ggsave(ggplot(ang, aes(group, V3))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig2.pdf')
ggsave(ggplot(ang, aes(group, V4))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig3.pdf')
ggsave(ggplot(ang, aes(group, V5))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig4.pdf')
ggsave(ggplot(ang, aes(group, V6))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig5.pdf')
ggsave(ggplot(ang, aes(group, V7))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig6.pdf')
ggsave(ggplot(ang, aes(group, V8))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig7.pdf')
ggsave(ggplot(ang, aes(group, V9))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig8.pdf')
ggsave(ggplot(ang, aes(group, V10))+geom_boxplot()+facet_grid(~split)+theme(axis.text.x = element_text(angle=90)), 
       filename = '/home/anne/Documents/paper_second_attempt.tex/figures/0.5/eig9.pdf')

