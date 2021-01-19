require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(RColorBrewer)
require(cowplot)
require(optparse)
theme_set(theme_cowplot())

option_list = list(
  make_option(c("-r", "--resultfolder"), action="store", default=NA, type='character',
              help="annoation file in gtf format"),
  make_option(c("-p", "--plotdir"), action="store", default=NA, type='character',
              help="The directory for plot"),
    make_option(c("-s", "--study.sizes"), action="store", default=NA, type='character',
              help="The directory for plot")
)
opt = parse_args(OptionParser(option_list=option_list))
plotdir <-opt$plotdir
setwd(opt$resultfolder)
alist<-list()

study_sizes<-fread(opt$study.sizes, header = T)

print(head(study_sizes))
for(set in list.dirs(recursive = F)){
  if(!(basename(set) %in% study_sizes$project.project_id)){
    next
  }
  a1<-fread(file.path(set, '0.5','single_site_bor_single_site_deflation_angles.tsv'), header=F)
  a1$comp <-'Eigendecompostion\n vs. Power Iteration (Deflation)'
  a2<-fread(file.path(set, '0.5','single_site_bor_single_site_subspace_angles.tsv'), header=F)
  a2$comp <-'Eigendecompostion\n vs. Subspace Iteration'
  a3<-fread(file.path(set, '0.5','single_site_subspace_single_site_deflation_angles.tsv'), header=F)
  a3$comp <-'Power Iteration (Deflation)\n vs. Subspace Iteration'
  a<-rbind(a1,a2,a3, fill=T)
  alist[[set]]<-a
}
all<-rbindlist(alist, fill=T)
print(colnames(all))
all<- all %>% pivot_longer(-c(V1, comp),names_to = 'rank', values_to = 'angle')
all<-as.data.table(all)
all<-all[!is.na(all$rank)]
all$rank <-as.factor(as.numeric(sapply(all$rank,function (x) gsub('V', "", x )))-1)
col<-c(brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 9),brewer.pal(name='Blues', n = 3))
g<-ggplot(all[angle<6], aes(y = angle, x = '', fill = rank))+geom_boxplot()+
  theme(plot.title = element_text(size = 35, hjust = 0.5), 
        legend.key.size = unit(0.7, "cm"),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 20), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=35),
        plot.subtitle = element_text(size = 25, hjust = 0.5),
        legend.text=element_text(size=20))+facet_wrap(~comp)+ylab('Angle [degree]')+
  scale_fill_manual('Rank', values = col )
g
ggsave(g, file = file.path(plotdir, 'angles_single_site_all_vs_all.eps'), dpi = 'print',
       width = 15, height = 10)
ggsave(g, file = file.path(plotdir, 'angles_single_site_all_vs_all.png'), dpi = 350,
       width = 15, height = 10)