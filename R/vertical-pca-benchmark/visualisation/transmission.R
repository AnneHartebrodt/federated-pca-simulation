require(data.table)
require(ggplot2)
require(tidyr)
require(cowplot)
require(optparse)
require(dplyr)
require(R.utils)
require(stringr)

palette_div <-
  c(
    '#062865',
    '#2f497d',
    '#4c6d96',
    '#6793af',
    '#81bac8',
    '#ffc4b3',
    '#f68888',
    '#d75161',
    '#ab203f',
    '#720022'
  )
palette_seq <-
  c(
    '#062865',
    '#203a72',
    '#324d80',
    '#43618d',
    '#52759b',
    '#618aa9',
    '#70a0b7',
    '#7eb6c5',
    '#8dcdd4',
    '#9ce4e2'
  )

my_theme <-
  theme_classic() + theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

byte2giga <- 100000000
byte2mega <- 100000


read_files_from_dir <- function(directory, prefix, suffix) {
  myfiles <-  list.files(directory)
  mytransmissions<- list()

  for (i in c(1,2, 3, 5, 10)) {
    for (f in which(str_detect(myfiles, paste0("^", prefix, i, ".*.", suffix)))) {
      filename <- myfiles[f]
      print(file.path(dir, filename))
      if(countLines(file.path(directory, filename))<=9){
        next
      }
      transmission <-
        fread(
          file.path(dir, filename),
          sep = '\t',
          skip = 9,
          header = F
        )
      transmission$file<-filename
      mytransmissions[[filename]]<-transmission
    }
  }
  mytransmissions<-rbindlist(mytransmissions)
  colnames(mytransmissions)<-c('action', 'iteration', 'client','rank', 'size', 'file')
  return (mytransmissions)
}
avg_tot_trans_per_cli<-function(transmissions){
  transmissions<-transmissions[client==1]
  avg_iterations<-ceiling(mean(transmissions[,max(iteration), by=.(file,rank)]$V1))
  print(avg_iterations)
  average_transmission<-round(sum(transmissions[,  mean(size), by=.(action)]$V1)/byte2giga, 4)
  print(average_transmission)
  total <- avg_iterations * average_transmission
  return(c(total, avg_iterations, average_transmission))
}
avg_vettor<-function(transmissions){
  transmissions<-transmissions[client==1 & rank<=10]
  #"G_i=SC"             "H_local=CS"         "H_global=SC"        "local_norm=CS"      "ALT_G_i_local=CS"  
  # "ALT_G_i=SC"         "global_norm=SC"     "local_dot_prod=CS"  "global_dot_prod=SC"
  transmissions<-transmissions[!(action %in% c('ALT_G_i_local=CS','ALT_G_i=SC'))]
  #ransmissions<-transmissions[action %in% c('ALT_G_i_local=CS','ALT_G_i=SC', "H_global=SC" , "H_local=CS")]
  # calculate the average # iterations per eigenvector
  # then take the sum of the average
  avg_iterations<-transmissions[,  max(iteration), by=.(file,rank)]
  avg_iterations<-avg_iterations[, mean(V1), by=.(rank)]
  avg_iterations<-ceiling(sum(avg_iterations$V1))
  print(avg_iterations)
  average_transmission<-round(sum(transmissions[,  mean(size), by=.(action)]$V1)/byte2giga, 4)
  print(average_transmission)
  total<-average_transmission*avg_iterations
  return(c(total, avg_iterations, average_transmission))
}


make_summary<-function(dir, prefix, suffix){
  print(paste0(prefix, '_'))
  transmissions<-read_files_from_dir(dir, paste0(prefix, suffix), 'transmission')

  tpc<-avg_tot_trans_per_cli(transmissions)
  print(tpc)
  
  transmissions<-read_files_from_dir(dir, paste0(prefix, '_guo', suffix), 'transmission')
  tpc_guo<-avg_vettor(transmissions)
  
  
  transmissions<-read_files_from_dir(dir, paste0(prefix, '_hybrid', suffix), 'transmission')
  tpc_hy<-avg_tot_trans_per_cli(transmissions)
  
  transmissions<-read_files_from_dir(dir, paste0(prefix,'_hybrid_reit', suffix), 'transmission')
  tpc_guo_hy<-avg_vettor(transmissions)
  
  total_hy <- tpc_hy+tpc_guo_hy
  
  df.sum<-as.data.table(t(data.table(tpc, tpc_guo, total_hy)))
  colnames(df.sum)<-c('Avg. transmission cost', 'Avg. Iterations', 'Cost per Iteration')
  df.sum$Algorithm<-c('Matrix', 'Vector', 'Hybrid')
  df.sum<-df.sum[, c(4,2,3,1)]
  return(df.sum)
}

# Mfeat data set
dir <-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mfeat/' 
prefix = 'mfeat\\-zer'
suffix = '_central_qr_'

df.sum.cent<-make_summary(dir, prefix, suffix)
df.sum.cent$qr<-'central'
df.sum.cent

suffix = '_fed_qr_'

df.sum<-make_summary(dir, prefix, suffix)
df.sum$qr<-'federated'

df.sum<-rbind(df.sum, df.sum.cent)
df.sum
fwrite(df.sum, 
       file = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/summaries/mfeat_transmission_summary.txt', 
       col.names = T, sep='&')

df.sum


# Mnist data set
dir <-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist/' 
prefix = 'raw'
suffix = '_central_qr_'
df.sum.cent<-make_summary(dir, prefix, suffix)
df.sum.cent$qr<-'central'
df.sum.cent

suffix = '_fed_qr_'
df.sum<-make_summary(dir, prefix, suffix)
df.sum$qr<-'federated'

df.sum<-rbind(df.sum, df.sum.cent)
df.sum
fwrite(df.sum, 
       file = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/summaries/mnist_transmission_summary.txt', 
       col.names = T, sep='&')



# Mnist data set
dir <-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS/' 
prefix = 'MMRF\\-COMMPASS'
suffix = '_central_qr_'
df.sum.cent<-make_summary(dir, prefix, suffix)
df.sum.cent$qr<-'central'
df.sum.cent

suffix = '_fed_qr_'
df.sum<-make_summary(dir, prefix, suffix)
df.sum$qr<-'federated'

df.sum<-rbind(df.sum, df.sum.cent)
df.sum
fwrite(df.sum, 
       file = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/summaries/MMRF-COMMPASS_transmission_summary.txt', 
       col.names = F, sep = '&')


# Mnist data set
dir <-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/1000g/chr1' 
prefix = 'chr1'
suffix = '_central_qr_'
df.sum.cent<-make_summary(dir, prefix, suffix)
df.sum.cent$qr<-'central'
df.sum.cent

suffix = '_fed_qr_'
df.sum<-make_summary(dir, prefix, suffix)
df.sum$qr<-'federated'

df.sum<-rbind(df.sum, df.sum.cent)
df.sum
fwrite(df.sum, 
       file = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/summaries/chr1_transmission_summary.txt', 
       col.names = T, sep='&')

df.sum

