# Run the PCA, outllier removal, PCA again
source('./code/federated_dp_pca/R/PCA_bor_aor_plot.R')


#Extract clusters with k-means
clusters<-kmeans(pca$x[,1:2], 2)
p.out<-ggbiplot(pca.outlier.free, var.axes = F, labels = lab, groups = as.factor(clusters$cluster))
p.out

cluster_result<-data.table(sample_id = 1:nrow(outlier.free), cluster_id = clusters$cluster)
