require(data.table)
require(ggplot2)
require(dplyr)


summarize_file<-function(path, nsites){
  if (file.exists(file.path(path, 'test_results/summary.tsv'))){
  summaries10<-fread(file.path(path, 'test_results/summary.tsv'), header=T, sep='\t')
  summaries10_2<-fread(file.path(path, 'run_summaries.tsv'), header=T, sep='\t')
  summaries10<-left_join(summaries10, summaries10_2)
  summaries10$nsites<-nsites
  head(summaries10)
  s10<-summaries10[ ,.(mean_runtime = round(mean(runtime),1), mean_iterations = round(mean(iterations),0), mean_bytes=round(mean(`total data sent`, na.rm=T),0) ,mean_angle_1 =round(mean(LSV1),2),mean_angle_5 =round(mean(LSV5),2),mean_angle_9 =round(mean(LSV9)),2), by=c('algorithm', 'nsites')]
  return(s10)
  }
  else{
    return NA
  }

}

s10<-summarize_file('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/app-tests/test-output/11/10/', 10)
s5<-summarize_file('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/app-tests/test-output/11/5/', 5)
s3<-summarize_file('/home/anne/Documents/featurecloud/pca/horizontal-pca/results/app-tests/test-output/11/3/', 3)

all<-rbind(s10,s5,s3)
all<- all[order(algorithm)]
fwrite(all,'/home/anne/Documents/featurecloud/pca/horizontal-pca/results/app-tests/test-output/11/stats.tsv', sep='\t')

mnist_summaries<-list()
counter <- 1
mnist_path <- '/home/anne/Documents/featurecloud/pca/horizontal-pca/results/app-tests/mnist-output/'
for (s in c(3,5,10)){
  for(i in 11:15){
    fp<- file.path(mnist_path, i, s, 'single')
    sum <-summarize_file(fp, s)
    if (!is.na(sum)){
    sum$seed<-i
    mnist_summaries[[counter]]<-sum
   counter<-counter+1
    }
  }
}


mnist_summaries<- rbindlist(mnist_summaries)
mnist_summaries<-mnist_summaries[ ,.(mean_runtime = round(mean(mean_runtime),1), mean_iterations = round(mean(mean_iterations),0), 
                     mean_bytes=round(mean(mean_bytes, na.rm=T),0) ,mean_angle_1 =round(mean(mean_angle_1),2),
                     mean_angle_5 =round(mean(mean_angle_5),2),mean_angle_9 =round(mean(mean_angle_9)),2), by=c('algorithm', 'nsites')]
mnist_summaries<-mnist_summaries[order(algorithm)]
fwrite(mnist_summaries,file.path(mnist_path, 'stats.tsv'), sep='\t')
