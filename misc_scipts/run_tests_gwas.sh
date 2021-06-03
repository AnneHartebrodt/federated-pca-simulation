gwaspath=$1
export PYTHONPATH=$PYTHONPATH:$gwaspath/federated_dp_pca
echo $PYTHONPATH
conda activate federated-pca
datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g

mkdir -p $resultpath
# take 2 chromosomes, we don't want to spam
for e in {1..2} ;
do
mkdir -p $resultpath/chr${e}
python3 $gwaspath/federated_dp_pca/python/PCA/vertical/vertical_pca_benchmark.py -f \
$datapath/chr${e}/chr${e}.thin \
--filetype 'gwas' --center -o $resultpath/chr${e} -r 10 -k 10 \
 -i 2000 --sep '\t' --header 0 --rownames 0 --names chr${e} --scale \
 --vert -s 5 --ortho_freq 100
done
echo "summaries"
for e in {1..2} :
do
cd $resultpath/chr${e}
echo $(pwd)
bash $gwaspath/federated_dp_pca/misc_scipts/make_summaries.sh $gwaspath/federated_dp_pca
outfile1='angles_precomp.tsv'
outfile2='angles_precomp.summary.tsv'
colname='angle'
Rscript $gwaspath/federated_dp_pca/R/vertical-pca-benchmark/data_cleanup/read_data.R -b . -s 'angles_precomp' -c $colname -o $outfile1 -d $gwaspath/federated_dp_pca
Rscript $gwaspath/federated_dp_pca/R/vertical-pca-benchmark/data_cleanup/aggregate_data.R -f $outfile1 -o $outfile2  -c $colname
Rscript $gwaspath/federated_dp_pca/R/vertical-pca-benchmark/data_cleanup/aggregate_data_with_dummy.R -f $outfile1 -o $outfile2  -c $colname
mkdir -p summaries
mv wide.* summaries
mv angles* co* eig* summaries
tar cvzf summaries.tar.gz summaries
done
