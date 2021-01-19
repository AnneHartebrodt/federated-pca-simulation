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



byte2giga <- 100000000
byte2mega <- 100000


read_files_from_dir <- function(directory, prefix, suffix) {
  myfiles <-  list.files(directory)
  myangles <- list()
  mylastlines <- list()
  for (i in c(1, 2, 3, 5, 10)) {
    for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix)))) {
      filename <- myfiles[f]
      print(file.path(directory, filename))
      if(countLines(file.path(directory, filename))<=9){
        next
      }
      angles <-
        fread(
          file.path(directory, filename),
          sep = '\t',
          skip = 9,
          header = F
        )
      if (is.null(angles)) {
        next
      }
      colnames(angles) <-
        c('iterations', sapply(1:(ncol(angles) - 1), function(x)
          paste0('', x)))
      angles$dataset <- filename
      angles$splits <- i
      
      angles <-
        as.data.table(pivot_longer(angles,-c(iterations, dataset, splits)))
      myangles[[filename]] <- angles
      
      
    }
  }

  myangles <- rbindlist(myangles)
  myangles <- myangles[!is.na(value)]
  
  return (myangles)
}




# dataset <- 'MMRF-COMMPASS'
# prefix = 'MMRF\\-COMMPASS_'
# s = 2
# offset <- c(rep(-50, 7),80 ,80)
# adjust <- 100

# dataset <- 'mnist'
# prefix = 'raw_'
# s = 2
# offset <- c(-15, -15, 0, -15,-15, -15,0 ,-15, 0)
# adjust <- 10


#dataset <- 'mfeat'
#prefix = 'mfeat\\-zer_'
#s = 1
#offset <- c(-5, -5, -5, 0,5, 0, -5 ,-5, -5)
#adjust = 5

dataset <- '1000g/chr1'
prefix = 'chr1_fed_qr_'
s = 2
offset <- rep(9, 0)
adjust <- 20

dirname <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/results/', dataset)

evals <-read_files_from_dir(dirname, prefix, 'eigenval')
outfile <- paste0('/home/anne/Documents/featurecloud/pca/vertical-pca/figures/', dataset,'_eigengaps_split_2.pdf')

evals[, diff:=abs(value-shift(value, 1)), by =.(iterations, dataset, splits)]
evals[, na:= paste0(shift(name, 1), '-' ,name), by =.(iterations, splits)]
evals$ldiff<-log(evals$diff)
mm <-mean(evals$ldiff, na.rm=T)
ssd <- sd(evals$ldiff, na.rm=T)

sum <- evals[,mean(diff), by =.(splits, name, na, iterations)]
sum<- sum[!is.na(V1)]


breaks <- sum %>% group_by(splits) %>%top_n(1, iterations) %>% select(na, splits, V1, iterations)
breaks<-as.data.table(breaks)
selection <-breaks$splits==s
breaks<-breaks[selection]

# hack the system here. Adjust the line label position
breaks$iterations<-breaks$iterations + offset
maxit <- breaks$iterations[1]-adjust


eigengaps <-ggplot(sum[splits==s & iterations<maxit], aes(iterations, V1, col=na))+
  geom_line()+
  scale_y_log10()+
  scale_color_manual(values=palette_div)+
  my_theme+
  guides(color=F)+
  ylab('Eigengap [Log10]')+
  geom_text(data = breaks, aes(label = na, x = iterations , y = V1, color = na))+
  theme(text = element_text(size=3))

eigengaps
ggsave(eigengaps, file=outfile, height = 10, units = 'cm')
