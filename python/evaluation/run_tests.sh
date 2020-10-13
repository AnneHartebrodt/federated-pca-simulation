#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist'
#export PYTHONPATH=$PYTHONPATH:/home/anne/Documents/featurecloud/pca/federated_dp_pca
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
#'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw' \
#--filetype 'mnist' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist' -r 20 -k 10 \
# -i 2000 --orthovector='/home/anne/Documents/featurecloud/pca/vertical-pca/results/mnist/orthogonal_vectors.angles'

#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mfeat'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
#'/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-fou,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-fac,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-kar,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-mor,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-pix,/home/anne/Documents/featurecloud/pca/vertical-pca/data/mfeat/mfeat-zer' \
#--filetype 'delim-list' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/mfeat' -r 20 -k 10 \
# -i 2000 --sep '\s+'


#mkdir -p '/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS'
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
#'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/output_transposed.txt' \
#--filetype 'delim' --center -o '/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS' -r 20 -k 10 \
# -i 2000 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale


gwaspath=$1
export PYTHONPATH=$PYTHONPATH:$gwaspath/federated_dp_pca
echo $PYTHONPATH
datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g
for e in {1..22} ; do
python3 $gwaspath/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
$datapath/chr${e}/chr${e}.thin  --scaled \
--filetype 'gwas' --center -o $resultpath/chr${e} -r 10 -k 10 \
 -i 2000 --sep '\t' --header 0 --rownames 0 --names chr${e} --variance --compare_pca $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.values
done

for e in {1..1} ; do
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
$datapath/chr${e}/chr${e}.thin  --scaled \
--filetype 'gwas' --center -o $resultpath/chr${e} -r 10 -k 10 \
 -i 2000 --sep '\t' --header 0 --rownames 0 --names chr${e} --variance --compare_pca $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.values
done
