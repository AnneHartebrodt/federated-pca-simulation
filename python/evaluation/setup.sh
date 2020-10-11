mkdir software
cd software
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20200923.zip
unzip plink2_linux_x86_64_20200923.zip
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip
unzip plink_linux_x86_64_20200616.zip

cd ..
git clone git@gitlab.com:hartebrodt/federated_dp_pca.git
git submodule update --init --recursive
conda env create -f federated_dp_pca/python/evaluation/environment.yml

#git clone git://github.com/samtools/htslib.git
#git clone git://github.com/samtools/bcftools.git
#cd bcftools
# The following is optional:
#autoheader && autoconf && ./configure --disable-libgsl --enable-perl-filters
#make
