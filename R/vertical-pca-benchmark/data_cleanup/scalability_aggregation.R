suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))

option_list = list(
  make_option(c("-b", "--basedir"), action="store", default=NA, type='character',
              help="base directory"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="output file"),
  make_option(c("-d", "--scriptdir"), action="store", default=NA, type='character',
              help="output file")
)
opt = parse_args(OptionParser(option_list=option_list))

scriptdir<-opt$scriptdir
source(file.path(scriptdir, "/R/vertical-pca-benchmark/data_cleanup/library.R"))


basedir<-opt$b
suffix <- opt$s
column_name <- opt$c
outfile<-opt$o
print(basedir)
#basedir<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/scalability2'

data_orientation<-c('0.1','0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9')
vector_or_subspace<-c('matrix', 'vector')
nr_sites<- c(2,3,5,10)
eigenvector_update<-c('power', 'gradient')
qr_method<-c('central_qr', 'federated_qr')


data <- read_all_files(basedir, 'transmission')
data<- rbindlist(data)

colnames(data)<- c("operation","iteration","client","eigenvector","nr_floats","data_fraction","matrix","sites","eigenvector_update","qr_method","filename")

summary <- data[client==1] %>% group_by(eigenvector, eigenvector_update, data_fraction, matrix, sites, qr_method, filename) %>% 
  summarise(max_it = max(iteration), sum_nr_floats = sum(nr_floats))

sumsum<-summary %>% group_by(eigenvector_update, data_fraction, matrix, sites, qr_method, filename) %>% summarise(max_it = sum(max_it), nr_float=sum(sum_nr_floats))
sumsum <- sumsum %>% group_by(eigenvector_update, data_fraction, matrix, sites, qr_method) %>% mutate(counter = row_number(matrix))
sumsum<-as.data.table(sumsum)

selection<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr","vector_3_gradient_central_qr")


wide.for.tikz<-sumsum %>% pivot_wider(id_cols = c(counter), names_from = c(matrix, sites, eigenvector_update, qr_method, data_fraction), values_from = c(max_it))
wide.for.tikz<-as.data.table(wide.for.tikz)
cols <- grep("matrix_3_power_central_qr|matrix_3_power_federated_qr|vector_3_gradient_central_qr", names(wide.for.tikz), value = TRUE)
cols<-c('counter', cols)
wide.for.tikz <- wide.for.tikz %>% select(cols)
fwrite(wide.for.tikz, paste0('wide.tikz.iterations.', outfile), sep='\t')

# make a summary for the maximal iteration
wide.maxit<-sumsum %>% pivot_wider(id_cols = c(data_fraction, counter), names_from = c(matrix, sites, eigenvector_update, qr_method), values_from = c(max_it))
wide.maxit<-as.data.table(wide.maxit)
long.maxit<-wide.maxit %>% pivot_longer(-data_fraction)
long.maxit<-as.data.table(long.maxit)
ggplot(long.maxit[name %in% selection], aes(as.factor(data_fraction), value, fill=name))+geom_boxplot()


wide.for.tikz<-sumsum %>% pivot_wider(id_cols = c(counter), names_from = c(matrix, sites, eigenvector_update, qr_method, data_fraction), values_from = c(nr_float))
wide.for.tikz<-as.data.table(wide.for.tikz)
cols <- grep("matrix_3_power_central_qr|matrix_3_power_federated_qr|vector_3_gradient_central_qr", names(wide.for.tikz), value = TRUE)
cols<-c('counter', cols)
wide.for.tikz <- wide.for.tikz %>% select(cols)
fwrite(wide.for.tikz, paste0('wide.tikz.transmission.', outfile), sep='\t')

# make a summary for the transmission cost
wide.transmission<-sumsum %>% pivot_wider(id_cols = c(counter), names_from = c(matrix, sites, eigenvector_update, qr_method, data_fraction), values_from = c(nr_float))
wide.transmission<-as.data.table(wide.transmission)
long.transmission<-wide.transmission %>% pivot_longer(-data_fraction)
long.transmission<-as.data.table(long.transmission)
ggplot(long.transmission[name %in% selection], aes(as.factor(data_fraction), value, fill=name))+geom_boxplot()+scale_y_log10()
