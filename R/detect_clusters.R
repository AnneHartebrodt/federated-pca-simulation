# Run the PCA, outllier removal, PCA again
source('./code/federated_dp_pca/R/PCA_bor_aor_plot.R')


#Extract clusters with k-means
clusters<-kmeans(pca.outlier.free$x[,1:2], 2)
p.out.clu<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab, groups = as.factor(clusters$cluster))

cluster_result<-data.table(sample_id = 1:nrow(outlier.free), cluster_id = clusters$cluster)

fwrite(cluster_result,
       file = file.path(opt$o, paste0(basename(dirname(opt$f)), '_clusters.tsv')), 
       col.names = F, row.names = F, quote = F, sep = '\t')

file.name.out<-paste0(basename(dirname(opt$f)),'_outlier_free_with_cluster', '.pdf')
ggsave(p.out.clu, filename = file.path(opt$o, file.name.out))
       