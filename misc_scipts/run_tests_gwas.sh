gwaspath=$1
export PYTHONPATH=$PYTHONPATH:$gwaspath/federated_dp_pca
echo $PYTHONPATH
conda activate federated-pca
datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g

# take 2 chromosomes, we don't want to spam
for e in {1..2} ;
do
python3 $gwaspath/federated_dp_pca/python/PCA/vertical/vertical_pca_benchmark.py-f \
$datapath/chr${e}/chr${e}.thin \
--filetype 'gwas' --center -o $resultpath/chr${e} -r 10 -k 10 \
 -i 200 --sep '\t' --header 0 --rownames 0 --names chr${e} --scale \
 --compare_pca $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.values \
 --vert -s 2,5
done
echo "summaries"
cd $resultpath/chr${e}
bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh \
$gwaspath/federated_dp_pca