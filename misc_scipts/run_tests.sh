basepath=$1
scriptpath=$2
export PYTHONPATH=$PYTHONPATH:$scriptpath

outfile=$basepath'/results-new-tests/mnist'


#mkdir -p $outfile
#python3 $scriptpath/python/PCA/vertical/vertical_pca_benchmark.py -f $basepath'/data/mnist/raw' \
#--filetype 'mnist' --center -o $outfile -r 20 -k 10 \
# -i 2000 -s 5 --vert --ortho_freq 100

echo "making summaries"
cd $outfile
echo $(pwd)
bash $scriptpath/misc_scipts/make_summaries.sh $scriptpath

##COMMPASS
#outfile='/home/anne/Documents/featurecloud/pca/vertical-pca/results/MMRF-COMMPASS'
#mkdir -p $outfile
#python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/python/PCA/vertical/vertical_pca_benchmark.py -f \
#'/home/anne/Documents/featurecloud/data/tcga/data_clean/MMRF-COMMPASS/coding_only.tsv' \
#--filetype 'delim' --center -o $outfile -r 20 -k 10 \
# -i 200 --sep '\t' --header 0 --rownames 0 --names 'MMRF-COMMPASS' --scale --vert -s 5,2 --hor
#echo "making summaries"
#cd $outfile
#bash /home/anne/Documents/featurecloud/pca/federated_dp_pca/misc_scipts/make_summaries.sh
