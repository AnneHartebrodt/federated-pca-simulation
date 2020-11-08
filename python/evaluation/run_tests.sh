export PYTHONPATH=$PYTHONPATH:/home/anne/Documents/featurecloud/pca/federated_dp_pca

mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist'
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
--filetype 'mnist' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist' -r 20 -k 10 \
 -i 2000 -s 10,5,3,2
#
mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/mnist'
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
--filetype 'mnist' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist' -r 20 -k 10 \
 -i 2000 -s 10,5,3,2 --unequal '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/splits.tsv'

mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mfeat'
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-fou,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-fac,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-kar,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-mor,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-pix,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-zer' \
--filetype 'delim-list' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mfeat' -r 20 -k 10 \
 -i 2000 --sep '\s+' -s 10,5,3,2

##
#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
#'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
#--filetype 'delim' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS' -r 20 -k 10 \
# -i 2000 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale

#
mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/MMRF-COMMPASS'
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
--filetype 'delim' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/MMRF-COMMPASS' -r 20 -k 10 \
 -i 2000 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale \
  --unequal '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/splits.tsv'