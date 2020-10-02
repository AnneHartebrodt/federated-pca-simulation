## Transform vcf format to bed format and filter for linkage disequilibrium
# use plink1 to do so
cd /work/Documents/gwas/data/1000g/raw
mkdir -p /work/Documents/gwas/results/1000g


for e in {1..22} ; do
  # make chromosome result folder
mkdir /work/Documents/gwas/results/1000g/chr${e}
cd chr${e}
mv ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz chr${e}

# nonfounders takes all the samples into account when performing pca.

/work/Documents/gwas/plink --vcf ALL.chr${e}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed -out chr${e};
/work/Documents/gwas/plink2 --bfile chr${e} -indep-pairwise 50 5 0.5 --out chr${e}.ld.indep;
/work/Documents/gwas/plink2 --bfile chr${e} --exclude chr${e}.ld.indep.out --thin-count 100000 --make-bed --out chr${e}.thin;
## Run plink to obtain eigenvectors
/work/Documents/gwas/plink2 --bfile chr${e}.thin --pca approx --nonfounders --maf 0.01 --out ../../../../results/1000g/chr${e}/chr${e}.thin
cd ..
done

