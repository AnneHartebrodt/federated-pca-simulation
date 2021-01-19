export PYTHONPATH=$PYTHONPATH:/home/anne/Documents/featurecloud/pca/federated_dp_pca

outfile='/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist_january'
mkdir -p $outfile
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
--filetype 'mnist' --center -o $outfile -r 1 -k 10 \
 -i 2 -s 5,2 --vert
echo "making summaries"
cd $outfile
bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh

##COMMPASS
outfile='/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS_jan'
mkdir -p $outfile
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical/vertical_pca_benchmark.py -f \
'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
--filetype 'delim' --center -o $outfile -r 1 -k 10 \
 -i 200 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale --vert -s 5,2
echo "making summaries"
cd $outfile
bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh
