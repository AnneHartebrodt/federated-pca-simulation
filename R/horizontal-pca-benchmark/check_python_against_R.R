require(data.table)
require(ggplot2)
require(viridis)

scale<-function(x, max.data, min.data){
  x<-(x-min.data)/(max.data-min.data)
  return(x)
}

raw.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/out.trunc.txt'
scaled.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/scaled.txt'
projections.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/out.text'

projections.datafile.multi.nonoise<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/nonoise_projections.txt'
eigenvectors.noinoise.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/nonoise_eigenvectors.tsv'
eigenvalues.nonoise.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/nonoise_eigenvalues.tsv'

projections.datafile.multi<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/projections.txt'
eigenvectors.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/eigenvectors.tsv'
eigenvalues.datafile<-'/home/anne/Documents/featurecloud/data/TCGA/htseq/eigenvalues.tsv'

raw<-fread(raw.datafile)
raw<-raw[,2:ncol(raw)]

# center and scale data
means<-colMeans(raw)
vaari<-apply(raw, 2, function(x) var(x))
centered2<-as.matrix(raw)
centered2<-apply(centered2, 1, function(x) (x-means)/vaari)
sum(centered2<-t(centered2))

# scale data between 0 and 1
max.data<-max(centered2)
min.data<-min(centered2)
scaled.data<-apply(centered2, c(1,2), function(x) scale(x, max.data,min.data))
scaled.data<-as.data.table(scaled.data)

pca.r.scaled<-prcomp(scaled.data, center = F)
ggplot(as.data.table(pca.r.scaled$x), aes(x=-PC1, y = PC2))+geom_point()


# data as produced by python app
py.scaled<-as.matrix(fread(scaled.datafile, header = T))
py.scaled<-py.scaled[,2:ncol(py.scaled)]
sum(as.matrix(scaled.data)-py.scaled)

pca.py.scaled<-prcomp(py.scaled, center = F)
ggplot(as.data.table(pca.py.scaled$x), aes(x=-PC1, y = PC2))+geom_point()

# projections as produced by standalone pca in python
projections1<-fread(projections.datafile)
ggplot(projections1, aes(x=V1, y=V2))+geom_point()

# projections distributed pca without noise
projections.dist.nonoise.py<-fread(projections.datafile.multi.nonoise)
ggplot(projections.dist.nonoise.py, aes(V1,V2, col=as.factor(V5)))+geom_point()+scale_color_viridis_d()

# projections with noise
projections.dist.py<-fread(projections.datafile.multi)
ggplot(projections.dist.py, aes(V1,V2, col=as.factor(V5)))+
  geom_point()+scale_color_viridis_d()+facet_grid('V7')

# eigenvalues without noise
eigenvalues.nonoise<-fread(eigenvalues.nonoise.datafile)
eigenvalues.r.scaled<-pca.r.scaled$sdev^2
# eigenvalues with noise
eigenvalues<-fread(eigenvalues.datafile)
eigenvalues.melted<-melt(eigenvalues, id= c('V5', 'V6', 'V7', 'V8'))
ggplot(eigenvalues.melted, aes(as.factor(V7), value))+geom_point()+
  geom_line()+facet_grid(cols = vars(V5), rows = vars(variable))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#eigenvectors
eigenvectors<-fread(eigenvectors.datafile)

