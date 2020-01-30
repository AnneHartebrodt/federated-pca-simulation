# Run the PCA, outllier removal, PCA again
source('PCA_bor_aor_plot.R')


#Extract clusters with k-means
clusters<-kmeans(pca$x[,1:2], 2)
