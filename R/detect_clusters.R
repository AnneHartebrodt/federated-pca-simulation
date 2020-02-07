# Run the PCA, outllier removal, PCA again
source('./code/federated_dp_pca/R/PCA_bor_aor_plot.R')



#Extract clusters with k-means
cluster_split<-function(x, opt, pca.outlier.free){
  k <- as.numeric(opt$number.clusters)
if(k>3){
clusters<-kmeans(pca.outlier.free$c[,1:k], centers=k, iter = 10000)
}else{
clusters<-kmeans(pca.outlier.free$x[,c(1,x)], centers = k, iter = 10000)
} 
p.out.clu<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab, groups = as.factor(clusters$cluster))

  cluster_result<-data.table(sample_id = 1:nrow(outlier.free), cluster_id = clusters$cluster)

  fwrite(cluster_result,
       file = file.path(opt$o, paste0(basename(dirname(opt$f)),'_', x, '_clusters.tsv')), 
       col.names = F, row.names = F, quote = F, sep = '\t')

  file.name.out<-paste0(basename(dirname(opt$f)),'_', x,'_outlier_free_with_cluster', '.pdf')
  ggsave(p.out.clu, filename = file.path(opt$o, file.name.out))
}

cluster_split(2, opt, pca.outlier.free)
cluster_split(3, opt, pca.outlier.free)

