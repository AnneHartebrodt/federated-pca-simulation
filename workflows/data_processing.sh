#!bin/bash
datadir=/home/anne/Documents/featurecloud/data
tcga_folder=$datadir/tcga

metadata_folder=$tcga_folder/metadata/
supplement_folder="/home/anne/Documents/featurecloud/pca/horizontal-pca/paper_supplement"

scripts="/home/anne/Documents/featurecloud/pca/federated_dp_pca"
preprocessing_scripts="/home/anne/Documents/featurecloud/code/data_preprocessing"

data_folder=$tcga_folder/data
study_folder=$tcga_folder/cancer_type


sample_provenance=$metadata_folder/sample_provenance.tsv
manifest=$metadata_folder/gdc_manifest.2020-04-17.txt
casefile="$metadata_folder/cases.2021-02-11.json"

#echo $metadata_folder
#python $preprocessing_scripts/organize_TCGA_metadata.py -f $casefile -o $sample_provenance
#Rscript $preprocessing_scripts/TCGA_studies_by_cancertype.R -m $manifest -s $supplement_folder -o $metadata_folder -p $sample_provenance
#ls $metadata_folder/cancer_type >> $metadata_folder/cancer_type/files.txt

#echo $metadata_folder
# $1 File with list of files
for file in $(cat $metadata_folder/cancer_type/files.txt)
do
	echo $file
	
	#Rscript /home/anne/Documents/featurecloud/code/data_preprocessing/organize_data.R \
	#-f $metadata_folder/cancer_type/$file -o 'output.txt' -b $data_folder -d $study_folder 
	#Rscript ../../../featurecloud-test/data_preprocessing/TCGA_rename_columns.R -f 'output.txt' -m $2 -o 'output_named.txt'
done

#mkdir -p $study_folder
for folder in $(ls $study_folder)
do
	echo $study_folder/$folder
	cd $study_folder/$folder
	echo $metadata_folder/samples"$folder".tsv
	#Rscript /home/anne/Documents/featurecloud/code/data_preprocessing/transpose_data.R -f output.txt -o output_transposed.txt
	#rm output_named.txt
	#rm output.txt
	#rm final_*
	#Rscript $preprocessing_scripts/codingGenesfromGTF.R  -f output_transposed.txt -c $datadir/genome/coding_genes.tsv -o coding_only.tsv  &> workflow.log
	Rscript $preprocessing_scripts/create_site_data_frames.R -f coding_only.tsv -m $metadata_folder/cancer_type/samples"$folder".tsv
	cd ..
done
