gwaspath=$1
export PYTHONPATH=$PYTHONPATH:$gwaspath/federated_dp_pca
echo $PYTHONPATH
conda activate federated-pca
datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g
for e in {1..22} ;
do
python3 $gwaspath/federated_dp_pca/python/PCA/guo_vertical_runner_benchmark_edition.py -f \
$datapath/chr${e}/chr${e}.thin \
--filetype 'gwas' --center -o $resultpath/chr${e} -r 10 -k 10 \
 -i 2000 --sep '\t' --header 0 --rownames 0 --names chr${e} --scale --compare_pca $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.values
done