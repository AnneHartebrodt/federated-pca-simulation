suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(optparse))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))

option_list = list(
  make_option(c("-f", "--infile"), action="store", default=NA, type='character',
              help="infile"),
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="outfile"),
  make_option(c("-c", "--column"), action="store", default=NA, type='character',
              help="columnname")
  
  
)
opt = parse_args(OptionParser(option_list=option_list))

infile<-opt$infile
outfile<-opt$outfile
column<-opt$column
print(infile)

#infile<-'/home/anne/Documents/featurecloud/pca/vertical-pca/results/1000g/chr2/summaries/angles.u.tsv'
#column<-'angle'

data<-fread(infile)
data[, 9]<-as.numeric(unlist(data[, 9]))

#split into matrix and vector
# matrix is as it should e
data.ma <- data[matrix=='matrix']
data.ve<-data[matrix=='vector']



# add continuous iteration counts to the vector tale
data.ve<-data[matrix=='vector'] %>% 
  group_by(orientation, matrix, sites, eigenvector_update, filename, qr_method, orthonormalisation_skip) %>% 
  mutate(iterations = row_number(matrix))
data.ve<-as.data.table(data.ve)

# find the starting iteration for each eigenvector related to the start of the complete run
minit<-data.ve %>% group_by(orientation, matrix, sites, eigenvector_update, filename, qr_method, orthonormalisation_skip,rank)%>%
  summarise(minit=min(iterations)-1)
minit<-as.data.table(minit)

# expand such that every iteration from the start of the global run
# to the start of the actual eigenvector gets 90 degrees
minit<-minit %>% uncount(minit) %>%  
  group_by(orientation, matrix, sites, eigenvector_update, filename, qr_method, orthonormalisation_skip, rank) %>% 
  mutate(iterations = row_number(matrix))
minit<-as.data.table(minit)
minit$angle<-90

# find max iteration to fill the dataframe to get correct average
maxit<-data.ve %>%
  group_by(orientation, matrix, sites, eigenvector_update, qr_method,orthonormalisation_skip,  rank)%>%
  summarise(maxit_all=max(iterations), angle=min(angle))
maxit1<-as.data.table(maxit)

maxit<-data.ve %>%
  group_by(orientation, matrix, sites, eigenvector_update, qr_method, orthonormalisation_skip, rank, filename)%>%
  summarise(maxit=max(iterations))
maxit<-as.data.table(maxit)

maxit<-merge(maxit, maxit1)
maxit[, diff:=maxit_all-maxit]
maxit<- maxit %>% uncount(diff) %>% group_by(orientation, matrix, sites, eigenvector_update, qr_method,orthonormalisation_skip, rank, filename) %>% mutate(iterations = row_number(filename)+ maxit)
maxit<-as.data.table(maxit)

maxit<-maxit %>% mutate(maxit_all=NULL, maxit=NULL)
maxit<-as.data.table(maxit)


# find max iteration to fill the dataframe to get correct average
maxit.ma<-data.ma %>%
  group_by(orientation, matrix, sites, eigenvector_update, qr_method,orthonormalisation_skip, rank)%>%
  summarise(maxit_all=max(iterations), angle=min(angle))
maxit1<-as.data.table(maxit.ma)

maxit.ma<-data.ma %>%
  group_by(orientation, matrix, sites, eigenvector_update, qr_method, orthonormalisation_skip, rank, filename)%>%
  summarise(maxit=max(iterations))
maxit.ma<-as.data.table(maxit.ma)

maxit.ma<-merge(maxit.ma, maxit1)
maxit.ma[, diff:=maxit_all-maxit]
maxit.ma<- maxit.ma %>% uncount(diff) %>% group_by(orientation, matrix, sites, eigenvector_update, qr_method,orthonormalisation_skip, rank, filename) %>% mutate(iterations = row_number(filename)+ maxit)
maxit.ma<-as.data.table(maxit.ma)

maxit.ma<-maxit.ma %>% mutate(maxit_all=NULL, maxit=NULL)
maxit.ma<-as.data.table(maxit.ma)

# paste everything together again.
data<-rbind(data.ma, minit, data.ve, maxit, maxit.ma)

# summarise
summary<-data %>% 
  select(iterations,orientation, matrix, sites, eigenvector_update, filename, qr_method, orthonormalisation_skip,rank, angle) %>%
  group_by(iterations,orientation, matrix, sites, eigenvector_update, qr_method,orthonormalisation_skip, rank) %>%
  summarise(mean_value = mean(get(column)))
summary<-as.data.table(summary)

fwrite(summary, paste0('long.dummy.', outfile), sep='\t')
print('Done')
# 
# ggplot(summary[sites==5 & qr_method=='central_qr' & matrix=="matrix" & eigenvector_update=='power'], 
#        aes(iterations, mean_value, col=as.factor(rank)))+geom_line()
# 
# # create column wise output format
# wide.for.tikz <-summary[orientation=='vertical'] %>% pivot_wider(id_cols = c(iterations), names_from = c(matrix, sites, eigenvector_update, qr_method, rank), values_from = mean_value)
# wide.for.tikz<-as.data.table(wide.for.tikz)
# 
# cols <- grep("matrix_5_power_central_qr|matrix_5_power_federated_qr|vector_5_gradient_central_qr", names(wide.for.tikz), value = TRUE)
# cols<-c('iterations', cols)
# wide.for.tikz <- wide.for.tikz %>% select(cols)
# fwrite(wide.for.tikz, paste0('wide.tikz.', outfile), sep='\t')
# 
# #make wide to create names
# wide.vertical<-summary[orientation=='vertical'] %>% pivot_wider(id_cols = c(iterations, rank), names_from = c(matrix, sites, eigenvector_update, qr_method), values_from = mean_value)
# wide.vertical<-as.data.table(wide.vertical)
# 
# #make long for ggplot
# d <- wide.vertical %>% pivot_longer(-c(iterations, rank))
# d<-as.data.table(d)
# d<-d[!is.na(value)]
# 
# # select the correct configuration
# selection<-c("matrix_5_power_central_qr", "matrix_5_power_federated_qr","vector_5_gradient_central_qr")
# #selection3<-c("matrix_3_power_central_qr", "matrix_3_power_federated_qr", "vector_3_power_central_qr", "vector_3_power_federated_qr")

# make the plot
# angles.plot<-ggplot(d[name %in% selection], aes(iterations, value, col=as.factor(name)))+
#   geom_line()+facet_wrap(~rank, scales = 'free')+
#   my_theme+ylab('Mean angle [degree]')+
#   xlab('#Iterations')+
#   scale_color_manual('Configuration', values = palette_div)+
#   theme(axis.line=element_line(),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         legend.position = c(0, 0.5),
#         legend.justification = c("right", "top"),
#         legend.box.just = "right",
#         legend.margin = margin(0.25, 0.25, 0.25, 0.25),
#         legend.title = element_text(size=10))+
#   guides(color=guide_legend(keyheight = 0.5, title = element_text('Configuration', size = 8)))
# 
# angles.plot

                             
