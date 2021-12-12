gwaspath=$1
export PYTHONPATH=$PYTHONPATH:$gwaspath/federated_dp_pca
echo $PYTHONPATH
#conda activate federated-pca
datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g

mkdir -p $resultpath
# take 2 chromosomes, we don't want to spam
for e in {1..2} ;
do
for i in {1..3}:
do
mkdir -p $resultpath/chr${e}
java -jar scaling.jar $datapath/chr${e}/chr${e}.$i.thin.traw $datapath/chr${e}/chr${e}.$i.thin.traw.scaled 0
python3 $gwaspath/federated_dp_pca/python/PCA/vertical/approximate_vertical_pca_benchmark.py -f \
$datapath/chr${e}/chr${e}.$i.thin \
--filetype 'gwas' --center -o $resultpath/chr${e}.$i -r 1 -k 10 \
 -i 200 --sep '\t' --header 0 --rownames 0 --names chr${e}.$i \
 --vert -s 5 --ortho_freq 1000
done
done

