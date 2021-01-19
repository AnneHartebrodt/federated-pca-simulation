export PYTHONPATH=$PYTHONPATH:/home/anne/Documents/featurecloud/pca/federated_dp_pca

outfile='/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist_january'
mkdir -p $outfile
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
--filetype 'mnist' --center -o $outfile -r 2 -k 10 \
 -i 200 -s 10,5,3,2 --vert
cd $outfile
bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh
#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/mnist'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
#'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
#--filetype 'mnist' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist' -r 20 -k 10 \
# -i 2000 -s 10,5,3,2 --unequal '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/splits.tsv'

##
outfile='/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS_jan'
mkdir -p $outfile
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
--filetype 'delim' --center -o $outfile -r 2 -k 10 \
 -i 200 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale --vert
cd $outfile
bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh
#
#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/MMRF-COMMPASS'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
#'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
#--filetype 'delim' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/MMRF-COMMPASS' -r 20 -k 10 \
# -i 2000 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale \
#  --unequal '/home/anne/Documents/featurecloud/pca/vertical-pca/results/unequal/splits.tsv'


#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/seeds'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical_pca_benchmark.py -f \
#'/home/anne/Documents/featurecloud/dev/pca-tool/fed-pca-client/test_datasets/seeds/seeds_dataset.csv' \
#--filetype 'delim' -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/seeds' -r 1 -k 6 \
# -i 2000 --sep ',' --header 0 --names 'seeds'
