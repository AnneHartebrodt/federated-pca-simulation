## Download hapmap data
conda activate fed-gwas
mkdir /work/Documents/gwas/data/hapmap/raw
cd /work/Documents/gwas/data/hapmap/raw

wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2
bunzip2 -c hapmap3_r1_b36_fwd.qc.poly.tar.bz2 | tar xvf -

export PYTHONPATH=$PYTHONPATH:/work/Documents/federatedPCA/code
export PYTHONPATH=$PYTHONPATH:/work/Documents/gwas/federated-gwas-scripts
export PYTHONPATH=$PYTHONPATH:/work/Documents/featurecloud/pca/federated_dp_pca

cd /work/Documents/gwas/data/hapmap/raw
for pop in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI
do
mkdir $pop;
mv hapmap3_r1_b36_fwd.${pop}.qc.poly.recode.ped $pop
cd $pop;
## Transform ped/map format to bed/bim/fam
/work/Documents/gwas/plink --ped hapmap3_r1_b36_fwd.${pop}.qc.poly.recode.ped --map hapmap3_r1_b36_fwd.${pop}.qc.poly.recode.map --make-bed -out hapmap3_r1_b36_fwd.${pop}.qc.poly.recode;

## Filter for LD and downsample data and produce raw data wit SNPs in rows
/work/Documents/gwas/plink2 --bfile hapmap3_r1_b36_fwd.${pop}.qc.poly.recode -indep-pairwise 50 5 0.5 --out hapmap3_r1_b36_fwd.${pop}.ld.indep;
/work/Documents/gwas/plink2 --bfile hapmap3_r1_b36_fwd.${pop}.qc.poly.recode --exclude hapmap3_r1_b36_fwd.${pop}.ld.indep.prune.out --make-bed --out hapmap3_r1_b36_fwd.${pop}.qc.poly.recode.pruned;
/work/Documents/gwas/plink2 --bfile hapmap3_r1_b36_fwd.${pop}.qc.poly.recode.pruned --thin-count 100000 --make-bed --out hapmap3_r1_b36_fwd.${pop}.thin;
/work/Documents/gwas/plink2 --bfile hapmap3_r1_b36_fwd.${pop}.thin --recode A-transpose --out hapmap3_r1_b36_fwd.${pop}.thin;

# run plink2 PCA
mkdir -p /work/Documents/gwas/results/hapmap
/work/Documents/gwas/plink2 --bfile hapmap3_r1_b36_fwd.${pop}.thin --pca approx --nonfounders --maf 0.01 --out ../../../../results/hapmap/${pop}/hapmap3_r1_b36_fwd.${pop}.thin

#cut -d$'\t' -f1-6 --complement hapmap3_r1_b36_fwd.${pop}.thin.traw > hapmap3_r1_b36_fwd.${pop}.thin.traw.values
python3 /home/anne/Documents/featurecloud/pca/federated_dp_pca/import_export/gwas_import.py -f hapmap3_r1_b36_fwd.${pop}.thin.traw -o hapmap3_r1_b36_fwd.${pop}.thin.traw.scaled -m 0.01 -b hapmap3_r1_b36_fwd.${pop}.thin.bim
#python3 /work/Documents/federatedPCA/code/federated-pca-lib/python/PCA/guo_vertical_runner.py -f hapmap3_r1_b36_fwd.${pop}.thin.traw.scaled -o hapmap3_r1_b36_fwd.${pop}.thin.scaled -k 10 -s 2 -i 3000 -p /work/Documents/gwas/results/hapmap/$pop


# compare PCA versions
python3 /work/Documents/gwas/federated-gwas-scripts/scaling/compare_gwas_pca.py --plink hapmap3_r1_b36_fwd.${pop}.thin.eigenvec --qr guo_multi_site_eigenvector.tsv --scikit guo_single_site_eigenvector.tsv --outputfile hapmap3_r1_b36_fwd.${pop}.thin.angles.2

cd ..;
done






