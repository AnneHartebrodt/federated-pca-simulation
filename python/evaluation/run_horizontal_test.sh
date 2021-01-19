export PYTHONPATH=$PYTHONPATH:/home/anne/Documents/featurecloud/pca/federated_dp_pca
mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/horizontal-results/mnist'
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/horizontal_benchmark.py -f \
'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
--filetype 'mnist' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/horizontal-results/mnist' -r 20 -k 10 \
 -i 2000 -s 10,5,3,2