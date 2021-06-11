# Federated Principal Component Analysis

This repository contains additional information regarding federated principal component analysis. Library code for the simulation of federated PCA is available in [vertical](/python/PCA/vertical). The folder [misc_scripts](/misc_scipts) contains examples on how to run the code and a basic workflow. Dependencies are managed using a conda environment and can be installed using [environment.yaml](./environment.yaml) 

The specific use case for principal component analysis handeled by this version of principal compoenent analysis is patient stratification for medical data. This is a very special use case, where the patients are distributed over serveral sites, which cannot share their data. Unlike 'feature' reduction PCA, where one wants to compute a lower dimensional approximation of the feature-by-feature covariance matrix, the PCA in GWAS is computed on the sample-by-sample covariance matrix. The samples are distributed over several sites, therefore the covariance matrix cannot be computed directly. It is not desirable to compute this covariance matrix, because the number of samples in a GWAS cohort can become quite large. Therefore here, we describe a covariance free alternative.

[federated-pca](./federated_pca.png)

## Demonstration code can be found here
If you want to read code for the simulation of federated vertically partitioned PCA, here is a nice simulation script: [Federated PCA](../python/PCA/vertical/simulate_federated_vertically_partionned_pca.py)

There is also one for federated Gram Schmidt orthonormalisation: [Gram-Schmidt](../python/PCA/vertical/simulate_federated_qr_orthonormalisation.py)

## Citation
If you use this repository in your work please cite as: 

## A Web service is available
If you would like to try federated principal component analysis in practice, please visit https://federated.compbio.sdu.dk/. The web service allows you to set up a federated study, which includes using the web interface to set the parameters, and downloading a user client, which handles the communication on the client site. More detailed information can be found on the web site. Please note that it is experimental software.
