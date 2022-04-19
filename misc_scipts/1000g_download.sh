## Transform vcf format to bed format and filter for linkage disequilibrium
# use plink1 to do so

gwaspath=$1
plink1path=$2
plink2path=$3

datapath=$gwaspath/data/1000g/raw
resultpath=$gwaspath/results/1000g
mkdir -p $datapath
mkdir -p $resultpath

echo $resultpath
cd $datapath
for e in {1..2} ; do
# make chromosome result folder
#mkdir -p $resultpath/chr${e}/plink
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
#mkdir -p chr${e}
#mv ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz chr${e}
#mv ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi chr${e}
# you are in data chromosome now
cd chr${e}
#$plink1path --vcf ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out chr${e};
## LD pruning with plink default values
## change maf filtering to here.
#$plink2path --bfile chr${e} --maf 0.01 --rm-dup exclude-all --make-bed --out chr${e}.rmdup
#$plink2path --bfile chr${e}.rmdup --indep-pairwise 50 5 0.5 --out chr${e}.ld.indep;
## downsampling to 100,000 SNPs for reasonable runtime.
#
#if [ ! -e "chr${e}.ld.indep.out" ] ; then
#    touch "chr${e}.ld.indep.out"
#fi

#all
$plink1path --bfile chr${e}.rmdup --recode A-transpose --out chr${e}.thin
java -jar gwaspath/federated_dp_pca/import_export/scaling.jar $datapath/chr${e}/chr${e}.thin.traw $datapath/chr${e}/chr${e}.traw.scaled.all 6 0

tail -n 500000 $datapath/chr${e}/chr${e}.traw.scaled.all > $datapath/chr${e}/chr${e}.traw.scaled.500000
head -n 100000 $datapath/chr${e}/chr${e}.traw.scaled.all > $datapath/chr${e}/chr${e}.traw.scaled.100000

# cut the phenotypic information out.
#cut -d$'\t' -f1-6 --complement chr${e}.thin.traw > chr${e}.thin.traw.values

## Run plink to obtain eigenvectors
# nonfounders takes all the samples into account when performing pca
# we are not interested in biologically sensible results so --nonfounder should be fine
# maf 0.01 recommended
#$plink2path --bfile chr${e}.thin --pca approx --nonfounders --maf 0.01 --out $resultpath/chr${e}/plink/chr${e}.thin
#cut -d$'\t' -f1-2 --complement $resultpath/chr${e}/plink/chr${e}.thin.eigenvec > $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.values
#$plink2path --bfile chr${e}.thin --pca --nonfounders --maf 0.01 --out $resultpath/chr${e}/plink/chr${e}.thin.exact
#cut -d$'\t' -f1-2 --complement $resultpath/chr${e}/plink/chr${e}.thin.exact.eigenvec > $resultpath/chr${e}/plink/chr${e}.thin.eigenvec.exact.values

cd ..

